open Core_kernel.Std
open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std
open Commons

(* FUNCTIONS *)


(* to parse rna conf file *)

type rna_sample = {
  id : string ;
  species : string ;
  ref_species : string ;
  sample_fastq : string sample_fastq ;
  run_trinity : bool ;
  path_assembly : string ;
  run_apytram : bool ;
}

type config_rna_seq = rna_sample list

type configuration = {
  config_rna_seq : config_rna_seq ;
  apytram_samples: rna_sample list ;
  trinity_samples : rna_sample list ;
  all_ref_samples : rna_sample list ;
  families : string list;
  rna_conf_file : string ;
  species_tree_file : string ;
  alignments_dir : string ;
  seq2sp_dir : string ;
  outdir : string ;
  threads : int;
  memory : int;
}

let parse_fastq_path = function
  | "-" -> None
  | x -> Some x
  
let parse_orientation = function
  | "F" -> Left F
  | "R" -> Left R
  | "US" -> Left US
  | "RF" -> Right RF
  | "FR" -> Right FR
  | "UP" -> Right UP
  | _ -> failwith {| Syntax error in configuration file (orientation must be in ["F","F","RF","FR","US","UP"] |}

let parse_line_fields_of_rna_conf_file = function
  | [ id ; species ; ref_species ; path_fastq_single ; path_fastq_left ; path_fastq_right ; orientation ; run_trinity ; path_assembly ; run_apytram] ->
     let run_trinity = match run_trinity with
       | "yes" | "Yes" | "y" | "Y" -> true
       | "no" | "No" | "n" | "N" -> false
       | _ -> failwith {| Syntax error in configuration file (run_trinity must be "yes" or "no" |}
     in
     let run_apytram = match run_apytram with
       | "yes" | "Yes" | "y" | "Y" -> true
       | "no" | "No" | "n" | "N" -> false
       | _ -> failwith {| Syntax error in configuration file (run_apytram must be "yes" or "no" |}
     in

     let sample_fastq = match (parse_fastq_path path_fastq_single, 
                             parse_fastq_path path_fastq_left, 
                             parse_fastq_path path_fastq_right,
                             parse_orientation orientation) with
       | ( None, Some  _, Some _, Right o ) -> Paired_end (path_fastq_left, path_fastq_right, o)
       | ( Some  _ , None, None, Left o ) -> Single_end (path_fastq_single, o)
       | _ -> failwith (path_fastq_single ^ path_fastq_left ^ path_fastq_right ^ orientation)(*{| Syntax error in configuration file (path_fastq_single must be "-" if data are "paired-end" and path_fastq_left and path_fastq_right must be "-" if data are "single-end".|}*)
     in
     { id ; species ; ref_species ; sample_fastq ; run_trinity ; path_assembly ; run_apytram }
  | _ -> failwith "Syntax error in configuration file"

let parse_rna_conf_file path =
  In_channel.read_lines path
  |> List.tl_exn
  |> List.map ~f:(String.split ~on:'\t')
  |> List.map ~f:parse_line_fields_of_rna_conf_file



let families_of_alignments_dir alignments_dir =
  Sys.readdir alignments_dir
  |> Array.filter ~f:(fun f -> 
    if Filename.check_suffix f ".fa" || Filename.check_suffix f ".fasta" then
      true
    else
      (printf "Warning: %s is not a fasta file\n" f ; false)
    )
  |> Array.map ~f:(fun f -> fst (String.lsplit2_exn f ~on:'.')) (* Il y a obligatoirement un point dans le nom du fichier fasta *)
  |> Array.to_list 
    
    
let load_configuration rna_conf_file species_tree_file alignments_dir seq2sp_dir threads memory outdir = 
  let config_rna_seq = parse_rna_conf_file rna_conf_file in
  let filter_samples f = List.filter config_rna_seq ~f in
  let id_list = List.map config_rna_seq ~f:(fun s -> s.id) in
  if List.contains_dup id_list then
    failwith {|There are duplicate id in the first colum of the config file.|}
  else
    { 
      config_rna_seq;
      apytram_samples = filter_samples (fun s -> s.run_apytram);
      trinity_samples = filter_samples (fun s -> s.run_trinity);
      all_ref_samples = filter_samples (fun s -> s.run_apytram || s.run_trinity);
      families = families_of_alignments_dir alignments_dir;
      rna_conf_file ;
      species_tree_file ;
      alignments_dir ;
      seq2sp_dir ;
      threads = int_of_string threads ;
      memory = int_of_string memory ;
      outdir;
    }

module Amalgam = struct

type configuration_dir = [ `configuration ] directory

let parse_input { rna_conf_file ; species_tree_file ; alignments_dir ; seq2sp_dir } : configuration_dir workflow= 
       workflow ~np:1 ~version:9 [
            cmd "../bin/ParseInput.py"  [ string rna_conf_file ;
                                          string species_tree_file;
                                          string alignments_dir;
                                          string seq2sp_dir;
                                          ident dest ;
                                        ]
    ]

let ref_transcriptomes =
    selector ["R_Sp_transcriptomes"] 

let ref_seq_fam_links =
    selector ["R_Sp_Seq_Fam_links"]

let ref_fams species family =
    selector ["R_Sp_Gene_Families"; species ^ "." ^ family ^ ".fa"] 

let ali_species2seq_links family =
    selector ["Alignments_Species2Sequences" ; "alignments." ^  family ^ ".sp2seq.txt" ]


  
let norm_fasta { all_ref_samples ; threads ; memory } =
  List.map all_ref_samples ~f:(fun s ->
    let seq_type = "fq" in 
    let max_cov = 50 in 
    let fastq = sample_fastq_map Bistro.Workflow.input s.sample_fastq in
    (s, Trinity.read_normalization seq_type memory max_cov threads fastq)
  )


let trinity_assemblies_of_norm_fasta norm_fasta memory threads =
  List.filter_map norm_fasta ~f:(fun (s,norm_fasta) -> 
    let is_paired = sample_fastq_is_paired s.sample_fastq in
    if s.run_trinity then
      Some (s, Trinity.trinity_fasta ~is_paired ~full_cleanup:true ~memory ~threads norm_fasta)
    else
      None
    )

let parse_seqids = true
let dbtype = "nucl"

let blast_dbs_of_norm_fasta norm_fasta =
  List.filter_map norm_fasta ~f:(fun (s,norm_fasta) -> 
    if s.run_apytram then
      Some (s, BlastPlus.makeblastdb ~parse_seqids ~dbtype  s.species norm_fasta)
    else
      None
    )
    

let seq_dispatcher ?s2s_tab_by_family query query_species ref_transcriptome seq2fam : fasta workflow = 
       workflow ~version:6 [
            mkdir_p tmp;
            cmd "../bin/SeqDispatcher.py"  [ 
              option (flag string "--sp2seq_tab_out_by_family" ) s2s_tab_by_family;
              opt "-tmp" ident tmp ;
              opt "-q" dep query ;
              opt "-qs" string query_species ;
              opt "-t" dep ref_transcriptome ;
              opt "-t2f" dep seq2fam;
              opt "-out" seq [ dest ; string ("/Trinity." ^ query_species )] ;
            ]
    ]

let trinity_annotated_fams_of_trinity_assemblies configuration_dir   =
    List.map ~f:(fun (s,trinity_assembly) ->
      let r = 
        seq_dispatcher 
          ~s2s_tab_by_family:true
          trinity_assembly 
          s.species 
          (configuration_dir / ref_transcriptomes / selector [ s.ref_species ^ "_transcriptome.fa" ])
          (configuration_dir / ref_seq_fam_links / selector [ s.ref_species ^ "_Fam_Seq.fa" ])
      in
      (s, r)
  )
         

let parse_apytram_results apytram_annotated_ref_fams =
  let config = tmp // "config.tsv" in
  workflow ~version:3 [
    mkdir_p tmp;
    heredoc ~dest:config (
      List.map apytram_annotated_ref_fams ~f:(fun (s, f, w) ->
        seq ~sep:"\t" [ string s.species ; string f ; dep w ]
      )
      |> seq ~sep:"\n"
    ) ;
    cmd "../bin/Parse_apytram_results.py"  [ config ; dest ]
    ]
    


let seq_integrator
      ?realign_ali
      ?log
      ?(species_to_refine_list = []) 
      ~family
      ~trinity_fam_results_dirs
      ~apytram_results_dir
      ~alignment_sp2seq
      alignment
      : _ directory workflow
      =
     
      let get_trinity_file_list extension dirs =
        List.map  dirs ~f:(fun (s,dir) ->
            [ dep dir ; string ("/Trinity." ^ s.species ^ "." ^ family ^ "." ^ extension) ; string ","]
            )
        |> List.concat
        in
        
      let get_apytram_file_list extension dir =
          [ dep dir ; string ("/apytram." ^ family ^ "." ^ extension) ; string ","]
      in
        
      
      
      let trinity_fasta_list  =  get_trinity_file_list "fa" trinity_fam_results_dirs in
      let trinity_sp2seq_list  =  get_trinity_file_list "sp2seq.txt" trinity_fam_results_dirs in
      
      let apytram_fasta  =  get_apytram_file_list "fa" apytram_results_dir in
      let apytram_sp2seq  =  get_apytram_file_list "sp2seq.txt" apytram_results_dir in
      
      let sp2seq = List.concat [[dep alignment_sp2seq ; string "," ] ; trinity_sp2seq_list ; apytram_sp2seq ]  in
      let fasta = List.concat [trinity_fasta_list; apytram_fasta]  in

       workflow ~version:11 [
            cmd "../bin/SeqIntegrator.py"  [ 
              opt "-tmp" ident tmp;
              opt "-ali" string alignment ;
              opt "-fa" (seq ~sep:"") fasta;
              opt "-sp2seq" (seq ~sep:"") sp2seq  ; (* list de sp2seq delimited by comas *)
              opt "-out" seq [ dest ; string "/" ; string family] ;
              opt "-sptorefine" (seq ~sep:",") (List.map species_to_refine_list ~f:(fun sp -> string sp) );
            ]
    ]



let merged_families_of_families configuration configuration_dir trinity_annotated_fams apytram_results_dir =
  List.map configuration.families ~f:(fun family ->
    let trinity_fam_results_dirs=
      List.map configuration.trinity_samples ~f:(fun s ->
        (s , List.Assoc.find_exn trinity_annotated_fams s)
        )
        in

    let alignment = configuration.alignments_dir ^ "/" ^ family ^ ".fa"  in
    let alignment_sp2seq = configuration_dir / ali_species2seq_links family in
    let species_to_refine_list = List.map configuration.all_ref_samples ~f:(fun s -> s.species)  in

    (family, seq_integrator ~species_to_refine_list ~family ~trinity_fam_results_dirs ~apytram_results_dir ~alignment_sp2seq  alignment )
    )


let merged_families_distributor merged_families =
  let config = tmp // "config.tsv" in
  let extension_list = [(".fa","Merged_fasta");(".tree","Merged_tree");(".sp2seq.txt","Sp2Seq_link")] in
  workflow ~version:1 [
    mkdir_p tmp;
    mkdir_p (dest // "Merged_fasta");
    mkdir_p (dest // "Merged_tree");
    mkdir_p (dest // "Sp2Seq_link");
    heredoc ~dest:config (
      List.map extension_list ~f:(fun (ext,dir) ->
          List.map  merged_families ~f:(fun (f, w) ->
          let input = w / selector [ f ^ ext ] in
          let output = dest // dir // (f ^ ext)  in
            seq ~sep:" " [ string "ln -s"; dep input ; ident output ]
          )
      |> seq ~sep:"\n"
      )
    |> seq ~sep:"\n"
    ) ; 
    cmd "bash" [ config ] 
    ]


let phyldog_of_merged_families_dirs configuration merged_families_dirs =
    let seqdir = merged_families_dirs / selector [ "Merged_fasta" ] in
    let treedir = merged_families_dirs / selector [ "Merged_tree" ] in
    let linkdir = merged_families_dirs / selector [ "Sp2Seq_link" ] in 
    let treefile = configuration.species_tree_file in
    let threads = configuration.threads in 
    let memory = configuration.memory in
    Phyldog.phyldog ~threads ~memory ~topogene:true ~timelimit:9999999 ~treefile ~linkdir ~treedir seqdir

let main configuration =
    let configuration_dir = parse_input configuration in
    let norm_fasta = norm_fasta configuration in
    let trinity_assemblies = trinity_assemblies_of_norm_fasta norm_fasta configuration.memory configuration.threads in
    let trinity_annotated_fams = trinity_annotated_fams_of_trinity_assemblies configuration_dir trinity_assemblies in
    
    let blast_dbs = blast_dbs_of_norm_fasta norm_fasta in
    
    let apytram_annotated_ref_fams =
        let pairs = List.cartesian_product configuration.apytram_samples configuration.families in
        List.map pairs ~f:(fun (s, fam) ->
          let query = configuration_dir / ref_fams s.ref_species fam in
          let blast_db = List.Assoc.find_exn blast_dbs s in
          let db_type = sample_fastq_orientation s.sample_fastq in
          (s, fam, Apytram.apytram ~plot:true ~i:3 ~query db_type blast_db)
          ) in
          
        
    let apytram_results_dir = parse_apytram_results apytram_annotated_ref_fams in
       
    
    let merged_families = merged_families_of_families configuration configuration_dir trinity_annotated_fams apytram_results_dir in
    
    let merged_families_dirs = merged_families_distributor merged_families in
       
    let phyldog = phyldog_of_merged_families_dirs configuration merged_families_dirs in
    
    let open Bistro_app in 
    List.concat [
      [ [ configuration.outdir ; "configuration" ] %>  configuration_dir ] ;
      List.map norm_fasta ~f:(fun (s,norm_fasta) ->
        [ configuration.outdir ; "norm_fasta" ; s.id ^ "_" ^ s.species ^ ".fa" ] %> norm_fasta
        );
      List.map trinity_assemblies ~f:(fun (s,trinity_assemblies) ->
        [ configuration.outdir ; "trinity_assemblies" ; "Trinity_assemblies." ^ s.id ^ "_" ^ s.species ^ ".fa" ] %> trinity_assemblies
        );
      List.map trinity_annotated_fams ~f:(fun (s,trinity_annotated_fams) ->
        [ configuration.outdir ; "trinity_annotated_fams" ; s.id ^ "_" ^ s.species ^ ".vs." ^ s.ref_species  ] %> trinity_annotated_fams
        );
      List.map blast_dbs ~f:(fun (s,blast_db) ->
        [ configuration.outdir ; "blast_db" ; s.id ^ "_"^ s.species ] %> blast_db
        );
      List.map apytram_annotated_ref_fams ~f:(fun (s, fam, apytram_result) ->
        [ configuration.outdir ; "apytram_annotated_fams" ; fam ; s.id ^ "_" ^ s.species ] %> apytram_result
        );
     [ [configuration.outdir ; "apytram_results" ] %> apytram_results_dir] ;
     
     [ [configuration.outdir ; "merged_families_dir" ] %> merged_families_dirs] ;
     
     [ [configuration.outdir ; "phyldog" ] %> phyldog] ;
     
    ]
 

end

(* RUN *)

let rna_conf_file = Sys.argv.(1)
let species_tree_file = Sys.argv.(2)
let alignments_dir = Sys.argv.(3)
let seq2sp_dir =  Sys.argv.(4)
let threads =  Sys.argv.(5)
let memory =  Sys.argv.(6)
let outdir =  Sys.argv.(7)


let configuration = load_configuration rna_conf_file species_tree_file alignments_dir seq2sp_dir threads memory outdir

let target_amalgam = Amalgam.main configuration 

let _ = Bistro_app.local ~np:configuration.threads  ~mem:( 1024 * configuration.memory) target_amalgam ~outdir

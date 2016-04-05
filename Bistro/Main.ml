open Core_kernel.Std
open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std

(* FUNCTIONS *)


(* to parse rna conf file *)
type rna_sample = {
  id : string ;
  species : string ;
  ref_species : string ;
  path_fastq : string ;
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
}

let parse_line_fields_of_rna_conf_file = function
  | [ id ; species ; ref_species ; path_fastq ; run_trinity ; path_assembly ; run_apytram] ->
     let run_trinity = match run_trinity with
       | "yes" | "Yes" -> true
       | "no" | "No" -> false
       | _ -> failwith "Syntax error in configuration file" 
     in
     let run_apytram = match run_apytram with
       | "yes" | "Yes" -> true
       | "no" | "No" -> false
       | _ -> failwith "Syntax error in configuration file" 
     in
     { id ; species ; ref_species ; path_fastq ; run_trinity ; path_assembly ; run_apytram }
  | _ -> failwith "Syntax error in configuration file"

let parse_rna_conf_file path =
  In_channel.read_lines path
  |> List.tl_exn
  |> List.map ~f:(String.split ~on:'\t')
  |> List.map ~f:parse_line_fields_of_rna_conf_file



let families_of_alignments_dir dir =
  Sys.readdir dir
  |> Array.filter ~f:(fun f -> 
    if Filename.check_suffix f ".fa" || Filename.check_suffix f ".fasta" then
      true
    else
      (printf "pas content %s n'est pas un fasta \n" f ; false)
    )
  |> Array.map ~f:(fun f -> fst (String.lsplit2_exn f ~on:'.')) (* Il y a obligatoirement un point dans le nom du fichier fasta *)
  |> Array.to_list 
    
let load_configuration rna_conf_file species_tree_file alignments_dir seq2sp_dir = 
  let config_rna_seq = parse_rna_conf_file rna_conf_file in
  let filter_samples f = List.filter config_rna_seq ~f in
  { 
    config_rna_seq;
    apytram_samples= filter_samples (fun s -> s.run_apytram) ;
    trinity_samples = filter_samples (fun s -> s.run_trinity);
    all_ref_samples = filter_samples (fun s -> s.run_apytram || s.run_trinity);
    families = families_of_alignments_dir alignments_dir;
    rna_conf_file ;
    species_tree_file ;
    alignments_dir ;
    seq2sp_dir ;
  }
  
let apytram_species { apytram_samples } =
  List.map apytram_samples ~f:(fun s -> s.ref_species)
  |> List.dedup
  

(* to run normalization *)
let memory = 1
let max_cov = 40
let nb_cpu = 2
let seq_type = "fq"

let targets_normalization samples =
  let open Bistro_app in
  List.map samples ~f:(fun s ->
    let fastq = Bistro.Workflow.input s.path_fastq in
    let norm_fasta = Trinity.read_normalization seq_type memory max_cov nb_cpu fastq in
    [
      [ "tmp_amalgam" ; "norm_fasta"; s.species ^ ".norm.fa" ] %> norm_fasta  ;
    ]
  )
  










module Amalgam = struct

type configuration_dir = [ `configuration ] directory

let parse_input { rna_conf_file ; species_tree_file ; alignments_dir ; seq2sp_dir } : configuration_dir workflow= 
       workflow ~version:4 [
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

let ref_fams =
    selector ["R_Sp_Fams"]

let ali_seq_2_species_links =
    selector ["Alignments_Sequences2Species"]

  
let norm_fasta { all_ref_samples } =
  List.map all_ref_samples ~f:(fun s ->
    let fastq = Bistro.Workflow.input s.path_fastq in
    (s, Trinity.read_normalization seq_type memory max_cov nb_cpu fastq)
  )


let trinity_assemblies_of_norm_fasta norm_fasta =
  List.filter_map norm_fasta ~f:(fun (s,norm_fasta) -> 
    if s.run_trinity then
      Some (s, Trinity.trinity ~full_cleanup:true ~memory norm_fasta)
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
       workflow [
            cmd "../bin/SeqDispatcher.py"  [ 
              option (flag string "--s2s_tab_out_by_family" ) s2s_tab_by_family;
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
         


let main configuration =
    let configuration_dir = parse_input configuration in
    let norm_fasta = norm_fasta configuration in
    let trinity_assemblies = trinity_assemblies_of_norm_fasta norm_fasta in
    let trinity_annotated_fams = trinity_annotated_fams_of_trinity_assemblies configuration_dir trinity_assemblies in
    
    let blast_dbs = blast_dbs_of_norm_fasta norm_fasta in
       
     
    
    let open Bistro_app in 
    List.concat [
      [ [ "tmp_amalgam" ; "configuration" ] %>  configuration_dir ] ;
      List.map norm_fasta ~f:(fun (s,norm_fasta) ->
        [ "tmp_amalgam" ; "norm_fasta" ; s.species ^ ".fa" ] %> norm_fasta
        );
      List.map trinity_assemblies ~f:(fun (s,trinity_assemblies) ->
        [ "tmp_amalgam" ; "trinity_assemblies" ; "Trinity_assemblies." ^ s.species ^ ".fa" ] %> trinity_assemblies
        );
      List.map trinity_annotated_fams ~f:(fun (s,trinity_annotated_fams) ->
        [ "tmp_amalgam" ; "trinity_annotated_fams" ; s.species ^ ".vs." ^ s.ref_species  ] %> trinity_annotated_fams
        );
      List.map blast_dbs ~f:(fun (s,blast_db) ->
        [ "tmp_amalgam" ; "blast_db" ; s.species ] %> blast_db
        );
        
    ]
 

end

(* RUN *)

let rna_conf_file = Sys.argv.(1)
let species_tree_file = Sys.argv.(2)
let alignments_dir = Sys.argv.(3)
let seq2sp_dir =  Sys.argv.(4)


let configuration = load_configuration rna_conf_file species_tree_file alignments_dir seq2sp_dir

let  target_amalgam = Amalgam.main configuration 

let _ = Bistro_app.local target_amalgam

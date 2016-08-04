(*
# File: amalgam.ml
# Created by: Carine Rey
# Created on: March 2016
#
#
# Copyright 2016 Carine Rey
# This software is a computer program whose purpose is to assembly
# sequences from RNA-Seq data (paired-end or single-end) using one or
# more reference homologous sequences.
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
*)

open Core.Std
open Bistro.Std
open Bistro.EDSL
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
  run_transdecoder : bool ;
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
     let run_transdecoder = true in
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
       | ( None, Some  _, Some _, Right o ) -> Fastq_Paired_end (path_fastq_left, path_fastq_right, o)
       | ( Some  _ , None, None, Left o ) -> Fastq_Single_end (path_fastq_single, o)
       | _ -> failwith (path_fastq_single ^ path_fastq_left ^ path_fastq_right ^ orientation)(*{| Syntax error in configuration file (path_fastq_single must be "-" if data are "paired-end" and path_fastq_left and path_fastq_right must be "-" if data are "single-end".|}*)
     in
     { id ;
       species ;
       ref_species ;
       sample_fastq ;
       run_trinity ;
       run_transdecoder ;
       path_assembly ;
       run_apytram
     }
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
  else if Filename.is_relative species_tree_file then
    failwith {|amalgam needs the absolute path of the species tree.|}
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
      threads ;
      memory ;
      outdir;
    }

module Amalgam = struct

type configuration_dir = [ `configuration ] directory

let parse_input { rna_conf_file ; species_tree_file ; alignments_dir ; seq2sp_dir } : configuration_dir workflow=
       workflow ~np:1 ~version:9 [
            cmd "ParseInput.py"  [ string rna_conf_file ;
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


let fastq_to_fasta_convertion {all_ref_samples} configuration_dir =
    List.map all_ref_samples ~f:(fun s ->
    let sample_fastq = sample_fastq_map input s.sample_fastq in
    let sample_fastq_to_sample_fasta = function
        | Fastq_Single_end (w, o ) -> Fasta_Single_end ( Trinity.fastool w , o )
        | Fastq_Paired_end (lw, rw , o) -> Fasta_Paired_end ( Trinity.fastool lw, Trinity.fastool rw , o)
    in

    let sample_fasta = sample_fastq_to_sample_fasta sample_fastq in
    (s,sample_fasta)
  )

let normalize_fasta fasta_reads {threads ; memory } =
  List.map fasta_reads ~f:(fun (s,fasta_sample) ->
    let max_cov = 50 in
    let normalization_dir = Trinity.fasta_read_normalization max_cov ~threads ~memory fasta_sample in

    let norm_fasta_sample_to_normalization_dir normalization_dir = function
        | Fasta_Single_end (w, o ) -> Fasta_Single_end ( normalization_dir / selector ["single.norm.fa"] , o )
        | Fasta_Paired_end (lw, rw , o) -> Fasta_Paired_end ( normalization_dir / selector ["left.norm.fa"] , normalization_dir / selector ["right.norm.fa"], o )
    in
    (s, norm_fasta_sample_to_normalization_dir normalization_dir fasta_sample )
  )


let trinity_assemblies_of_norm_fasta norm_fasta { memory ; threads } =
  List.filter_map norm_fasta ~f:(fun (s, norm_fasta) ->
    if s.run_trinity then
      Some (s, Trinity.trinity_fasta ~full_cleanup:true ~memory ~threads norm_fasta)
    else
      None
    )

let transdecoder_orfs_of_trinity_assemblies trinity_assemblies { memory ; threads } =
  List.map trinity_assemblies ~f:(fun (s,trinity_assembly) ->
  if s.run_transdecoder  then
      let pep_min_length = 50 in
      let retain_long_orfs = 150 in
      (s, Transdecoder.transdecoder ~retain_long_orfs ~pep_min_length ~only_best_orf:true ~memory ~threads trinity_assembly)
    else
      (s, trinity_assembly)
    )

let concat = function
  | [] -> raise (Invalid_argument "fastX concat: empty list")
  | x :: [] -> x
  | fXs ->
    workflow ~descr:"fastX.concat" [
      cmd "cat" ~stdout:dest [ list dep ~sep:" " fXs ]
    ]

let parse_seqids = true
let dbtype = "nucl"

let blast_dbs_of_norm_fasta norm_fasta =
  List.filter_map norm_fasta ~f:(fun (s, norm_fasta) ->
    if s.run_apytram then
      let fasta_to_norm_fasta_sample = function
        | Fasta_Single_end (w, _ ) -> w
        | Fasta_Paired_end (lw, rw , _) -> concat [ lw ; rw ]
      in
      let fasta = fasta_to_norm_fasta_sample norm_fasta in
      Some (s, BlastPlus.makeblastdb ~parse_seqids ~dbtype  (s.id ^ "_" ^ s.species) fasta)
    else
      None
    )


let seq_dispatcher ?s2s_tab_by_family query query_species query_id ref_transcriptome seq2fam : fasta workflow =
       workflow ~version:6 [
            mkdir_p tmp;
            cmd "SeqDispatcher.py"  [
              option (flag string "--sp2seq_tab_out_by_family" ) s2s_tab_by_family;
              opt "-tmp" ident tmp ;
              opt "-q" dep query ;
              opt "-qs" string query_species ;
              opt "-qid" string query_id ;
              opt "-t" dep ref_transcriptome ;
              opt "-t2f" dep seq2fam;
              opt "-out" seq [ dest ; string ("/Trinity." ^ query_id ^ "." ^ query_species )] ;
            ]
    ]

let trinity_annotated_fams_of_trinity_assemblies configuration_dir   =
    List.map ~f:(fun (s,trinity_assembly) ->
      let r =
        seq_dispatcher
          ~s2s_tab_by_family:true
          trinity_assembly
          s.species
          s.id
          (configuration_dir / ref_transcriptomes / selector [ s.ref_species ^ "_transcriptome.fa" ])
          (configuration_dir / ref_seq_fam_links / selector [ s.ref_species ^ "_Fam_Seq.fa" ])
      in
      (s, r)
  )

let apytram_orfs_ref_fams_of_apytram_annotated_ref_fams apytram_annotated_ref_fams memory =
    List.map apytram_annotated_ref_fams ~f:(fun (s, f, apytram_result_fasta) ->
      if s.run_transdecoder then
        let pep_min_length = 20 in
        let retain_long_orfs = 150 in
        let filtered_orf = Transdecoder.transdecoder ~retain_long_orfs ~pep_min_length ~only_best_orf:true ~threads:1 ~memory apytram_result_fasta in
        (s, f, filtered_orf)
      else
        (s, f, apytram_result_fasta)
    )


let parse_apytram_results apytram_annotated_ref_fams =
  let config = Bistro.Expr.(
      List.map apytram_annotated_ref_fams ~f:(fun (s, f, w) ->
        seq ~sep:"\t" [ string s.species ; string s.id ; string f ; dep w ]
      )
      |> seq ~sep:"\n"
    )
  in
  workflow ~version:4 [
    script "Parse_apytram_results.py" config ~args:[dest]
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

      let tmp_merge = dest // "tmp" in

       workflow ~version:11 [
            mkdir_p tmp_merge ;
            cmd "SeqIntegrator.py"  [
              opt "-tmp" ident tmp_merge;
              opt "-ali" string alignment ;
              opt "-fa" (seq ~sep:"") fasta;
              option (flag string "--realign_ali") realign_ali;
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

    (family, seq_integrator ~realign_ali:true ~species_to_refine_list ~family ~trinity_fam_results_dirs ~apytram_results_dir ~alignment_sp2seq  alignment )
    )


let merged_families_distributor merged_families =
  let extension_list = [(".fa","Merged_fasta");(".tree","Merged_tree");(".sp2seq.txt","Sp2Seq_link")] in
  workflow ~version:1 [
    mkdir_p tmp;
    mkdir_p (dest // "Merged_fasta");
    mkdir_p (dest // "Merged_tree");
    mkdir_p (dest // "Sp2Seq_link");
    let config = Bistro.Expr.(
        List.map extension_list ~f:(fun (ext,dir) ->
            List.map  merged_families ~f:(fun (f, w) ->
                let input = w / selector [ f ^ ext ] in
                let output = dest // dir // (f ^ ext)  in
                seq ~sep:" " [ string "ln -s"; dep input ; ident output ]
              )
            |> seq ~sep:"\n"
          )
        |> seq ~sep:"\n"
      )
    in
    script "bash" config
  ]


let phyldog_of_merged_families_dirs configuration merged_families_dirs =
    let seqdir = merged_families_dirs / selector [ "Merged_fasta" ] in
    let treedir = merged_families_dirs / selector [ "Merged_tree" ] in
    let linkdir = merged_families_dirs / selector [ "Sp2Seq_link" ] in
    let treefile = configuration.species_tree_file in
    let threads_max = (List.length configuration.families) + 1 in
    let threads = Pervasives.min threads_max configuration.threads in
    let memory = configuration.memory in
    Phyldog.phyldog ~threads ~memory ~topogene:true ~timelimit:9999999 ~treefile ~linkdir ~treedir seqdir


let output_of_phyldog phyldog merged_families families =
    workflow ~version:1 [
      mkdir_p (dest // "Alignments");
      mkdir_p (dest // "Sp2Seq_link");
      mkdir_p (dest // "Gene_trees");
      let extension_list = [(".fa","Alignments");(".sp2seq.txt","Sp2Seq_link")] in
      let config = Bistro.Expr.(
        seq ~sep:"\n" [
          List.map extension_list ~f:(fun (ext,dir) ->
                List.map  merged_families ~f:(fun (f, w) ->
                  let input = w / selector [ f ^ ext ] in
                  let output = dest // dir // (f ^ ext)  in
                  seq ~sep:" " [ string "ln -s"; dep input ; ident output ]
                  )
                  |> seq ~sep:"\n"
            )
            |> seq ~sep:"\n" ;
          let (ext,dir) = (".ReconciledTree","Gene_trees/") in
          List.map families ~f:(fun f ->
                 let input = phyldog / selector [ dir ^ f ^ ext ] in
                 let output = dest // dir // (f ^ ".tree")  in
                 seq ~sep:" " [ string "ln -s"; dep input ; ident output ]
                 )
                 |> seq ~sep:"\n";
        ]
      )
      in
    script "bash" config;
    ]

let main configuration =
    let divided_memory = Pervasives.(configuration.memory / configuration.threads) in

    let configuration_dir = parse_input configuration in

    let fasta_reads = fastq_to_fasta_convertion configuration configuration_dir in

    let norm_fasta = normalize_fasta fasta_reads configuration in

    let trinity_assemblies = trinity_assemblies_of_norm_fasta norm_fasta configuration in

    let trinity_orfs = transdecoder_orfs_of_trinity_assemblies trinity_assemblies configuration in

    let trinity_annotated_fams = trinity_annotated_fams_of_trinity_assemblies configuration_dir trinity_orfs in

    let blast_dbs = blast_dbs_of_norm_fasta norm_fasta in

    let apytram_annotated_ref_fams =
        let pairs = List.cartesian_product configuration.apytram_samples configuration.families in
        List.map pairs ~f:(fun (s, fam) ->
          let query = configuration_dir / ref_fams s.ref_species fam in
          let blast_db = List.Assoc.find_exn blast_dbs s in
          let db_type = sample_fastq_orientation s.sample_fastq in
          let w = Apytram.apytram ~no_best_file:true ~write_even_empty:true ~plot:false ~i:5 ~memory:divided_memory ~query db_type blast_db in
          let apytram_filename = "apytram." ^ s.ref_species ^ "." ^ fam ^ ".fasta" in
          (s, fam, w / selector [ apytram_filename ] )
          )
        in

    let apytram_orfs_ref_fams = apytram_orfs_ref_fams_of_apytram_annotated_ref_fams apytram_annotated_ref_fams divided_memory in

    let apytram_results_dir = parse_apytram_results apytram_orfs_ref_fams in

    let merged_families = merged_families_of_families configuration configuration_dir trinity_annotated_fams apytram_results_dir in

    let merged_families_dirs = merged_families_distributor merged_families in

    let phyldog = phyldog_of_merged_families_dirs configuration merged_families_dirs in

    let output = output_of_phyldog phyldog merged_families configuration.families in



    let open Bistro_app in

    let target_to_sample_fasta s d = function
          | Fasta_Single_end (w, _ ) -> [[ d ; s.id ^ "_" ^ s.species ^ ".fa" ] %> w ]
          | Fasta_Paired_end (lw, rw , _) -> [[ d ; s.id ^ "_" ^ s.species ^ ".left.fa" ] %> lw ; [ d ; s.id ^ "_" ^ s.species ^ ".right.fa" ] %> lw]
     in

    List.concat [
      [[ "configuration" ] %>  configuration_dir ]
        ;
      List.concat (List.map fasta_reads ~f:(fun (s,sample_fasta) -> target_to_sample_fasta s "raw_fasta" sample_fasta))
        ;
      List.concat (List.map norm_fasta ~f:(fun (s,norm_fasta) -> target_to_sample_fasta s "norm_fasta" norm_fasta))
        ;
      List.map trinity_assemblies ~f:(fun (s,trinity_assembly) ->
        [ "trinity_assemblies" ; "Trinity_assemblies." ^ s.id ^ "_" ^ s.species ^ ".fa" ] %> trinity_assembly
        )
        ;
      List.map trinity_orfs ~f:(fun (s,trinity_orf) ->
        [ "trinity_assemblies" ; "Transdecoder_cds." ^ s.id ^ "_" ^ s.species ^ ".fa" ] %> trinity_orf
        )
        ;
      List.map trinity_annotated_fams ~f:(fun (s,trinity_annotated_fams) ->
        [ "trinity_annotated_fams" ; s.id ^ "_" ^ s.species ^ ".vs." ^ s.ref_species ] %> trinity_annotated_fams
        )
        ;
      List.map blast_dbs ~f:(fun (s,blast_db) ->
        [ "blast_db" ; s.id ^ "_" ^ s.species ] %> blast_db
        )
        ;
      List.map apytram_annotated_ref_fams ~f:(fun (s, fam, apytram_result) ->
        [ "apytram_annotated_fams" ; fam ; s.id ^ "_" ^ s.species ] %> apytram_result
        )
        ;
      List.map apytram_orfs_ref_fams ~f:(fun (s, fam, apytram_result) ->
        [ "apytram_transdecoder_orfs" ; fam ; s.id ^ "_" ^ s.species ] %> apytram_result
        )
        ;
      [["apytram_results" ] %> apytram_results_dir]
        ;
      List.map merged_families ~f:(fun (fam, merged_family) ->
        [ "merged_families" ; fam  ] %> merged_family
        )
        ;
      [["merged_families_dir" ] %> merged_families_dirs]
        ;
      [["phyldog" ] %> phyldog]
        ;
      [[ "output" ] %> output ]
        ;
    ]


end

(* RUN *)

let main config_file outdir species_tree_file alignments_dir seq2sp_dir np memory () =
  let np = Option.value ~default:1 np in
  let memory = Option.value ~default:1 memory in
  let configuration = load_configuration config_file species_tree_file alignments_dir seq2sp_dir np memory outdir in
  let target_amalgam = Amalgam.main configuration in
  Bistro_app.local ~np:configuration.threads  ~mem:( 1024 * configuration.memory) target_amalgam ~outdir

let spec =
  let open Command.Spec in
  empty
  +> flag "--config"          (required file)   ~doc:"PATH Configuration file."
  +> flag "--outdir"          (required string) ~doc:"PATH Destination directory."
  +> flag "--species-tree"    (required file)   ~doc:"ABSOLUTE PATH Species tree file in nw format containing all species. Warning absolute path is required."
  +> flag "--alignment-dir"   (required string) ~doc:"PATH Directory containing all gene family alignments (Family_name.fa) in fasta format."
  +> flag "--seq2sp-dir"      (required string) ~doc:"PATH Directory containing all link files (Family_name.tsv). A line for each sequence and its species spaced by a tabulation."
  +> flag "--np"              (optional int)    ~doc:"INT Number of CPUs. (Default:1)"
  +> flag "--memory"          (optional int)    ~doc:"INT Number of GB of system memory to use.(Default:1)"

let command =
  Command.basic
    ~summary:"Amalgam"
    spec
    main

let () = Command.run ~version:"0.1" command

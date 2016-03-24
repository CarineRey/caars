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
  run_apytram : bool
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


(* to run Trinity on RNA samples *)
let targets_trinity samples =
  let open Bistro_app in
  List.map samples ~f:(fun s ->
    let fastq = Bistro.Workflow.input s.path_fastq in
    let assembly = Trinity.trinity fastq in
    [
      [ "output" ; s.path_assembly ] %> assembly  ;
    ]
  )



(* to run makeblastdb on RNA samples *)

let parse_seqids = "yes"
let dbtype = "nucl"

let targets_apytram samples =
  let open Bistro_app in
  List.map samples ~f:(fun s ->
    let fastq = Bistro.Workflow.input s.path_fastq in
    let fasta = Trinity.fastq2fasta fastq in
    let db_blast = BlastPlus.makeblastdb dbtype s.species fasta in
    [
      [ "tmp" ; "fasta"; s.species ^ ".fa" ] %> fasta  ;
      [ "tmp" ; "BlastDB"; s.species ^ "_DB" ] %> db_blast  ;
    ]
  )


(* to check path exits *)
let check_path_exists path =
  if not (Sys.file_exists path) then (
    let msg = "No such file or directory " ^ path in
    failwith msg
  )

(* RUN *)

(* Parsing rna conf file *)
let rna_conf_file = Sys.argv.(1)
let species_tree_file = Sys.argv.(2)
let ali_dir = Sys.argv.(3)
let seq2sp_dir =  Sys.argv.(4)

let parsed_rna_conf_file = parse_rna_conf_file rna_conf_file

(* Check if all paths exist*)
let () =
  List.iter parsed_rna_conf_file ~f:(fun sample ->
      check_path_exists sample.path_fastq
    )

(** Get each used reference species **)
let used_ref_species = parsed_rna_conf_file
  |> List.map ~f:(fun s -> s.ref_species )
  |> Caml.List.sort_uniq compare
  
(* Run ParseInput.py *)


let parse_input config_file species_tree_file ali_dir seq2sp_dir : fasta workflow= 
       workflow [
            cmd "../bin/ParseInput.py"  [ string config_file ;
                                          string species_tree_file;
                                          string ali_dir;
                                          string seq2sp_dir;
                                          ident dest ;
                                        ]
    ]


let temporary_transcriptome_dir  = parse_input rna_conf_file species_tree_file ali_dir seq2sp_dir

let target_parse_input = let open Bistro_app in
                            [[ "tmp";"test"] %> temporary_transcriptome_dir; ]
let _ = Bistro_app.local target_parse_input

(* Read Normalization of RNA sample*)
(*
let _ = parsed_rna_conf_file
  |> List.filter ~f:(fun s -> let res = false in 
                                  (if (s.run_apytram) then
                                   let res = true )
                                res)
*)

(* Run makeblastdb on RNA samples *)
let _ =
  parsed_rna_conf_file
  |> List.filter ~f:(fun s -> s.run_apytram )
  |> targets_apytram
  |> List.concat
  |> Bistro_app.local

(* Run Trinity on RNA samples *)
let _  =
  parsed_rna_conf_file
  |> List.filter ~f:(fun s -> s.run_trinity )
  |> targets_trinity
  |> List.concat
  |> Bistro_app.local





(* Run Seq_Dispatcher.py on Trinity assemblies *)
(** For each used reference species build a fasta with all its sequences **)

let seq_dispatcher query target seq2sp species : fasta workflow = 
       workflow [
            cmd "../bin/SeqDispatcher.py"  [ string "-q" ; string query ;
                                              string "-t" ; string target ;
                                              string "-t2f"; string seq2sp;
                                              string "-out"; seq [ dest ; string ("/Trinity_" ^ species )] ;
                                        ]
    ]

(*let _ = List.iter ~f:(printf "%s\n" ) used_ref_species*)


(*let _ = print_string fasta_paths *)
(*let _ = List.iter ~f:(printf "\n%s\n" ) all_ref_seqs_fasta*)
(*let _ = print_string all_ref_seqs_fasta*)



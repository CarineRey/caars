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


let fastq2fasta (fastq : _ fastq workflow ) : fasta workflow = 
       workflow [
            cmd "awk" ~stdout:dest [ string "NR%4==1||NR%4==2" ; dep fastq ; string "| tr \"@\" \">\""]
    ]


let parse_seqids = "yes"
let dbtype = "nucl"
(* to run makeblastdb on RNA samples *)
let targets_apytram samples =
  let open Bistro_app in
  List.map samples ~f:(fun s ->
    let fastq = Bistro.Workflow.input s.path_fastq in
    let fasta = fastq2fasta fastq in
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
  

(* Run Trinity on RNA samples *)
let _  =
  parsed_rna_conf_file
  |> List.filter ~f:(fun s -> s.run_trinity )
  |> targets_trinity
  |> List.concat
  |> Bistro_app.local

(* Run makeblastdb on RNA samples *)
let _ =
  parsed_rna_conf_file
  |> List.filter ~f:(fun s -> s.run_apytram )
  |> targets_apytram
  |> List.concat
  |> Bistro_app.local



(* Run Seq_Dispatcher.py on Trinity assemblies *)
(** For each used reference species build a fasta with all its sequences **)


(*let _ = List.iter ~f:(printf "%s\n" ) used_ref_species*)


(*let _ = print_string fasta_paths *)
(*let _ = List.iter ~f:(printf "\n%s\n" ) all_ref_seqs_fasta*)
(*let _ = print_string all_ref_seqs_fasta*)



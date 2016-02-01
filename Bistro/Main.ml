open Core_kernel.Std
open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std




(* FUNCTIONS *)


(* Parsing rna conf file *)
type rna_sample = {
  id : string ;
  species : string ;
  ref_species : string ;
  path_fastq : string ;
  run_trinity : bool ;
  path_assembly : string ;
}

let parse_line_fields_of_rna_conf_file = function
  | [ id ; species ; ref_species ; path_fastq ; run_trinity ; path_assembly ] ->
     let run_trinity = match run_trinity with
       | "yes" | "Yes" -> true
       | "no" | "No" -> false
       | _ -> failwith "Syntax error in configuration file" 
     in
     { id ; species ; ref_species ; path_fastq ; run_trinity ; path_assembly }
  | _ -> failwith "Syntax error in configuration file"

let parse_rna_conf_file path =
  In_channel.read_lines path
  |> List.tl_exn
  |> List.map ~f:(String.split ~on:'\t')
  |> List.map ~f:parse_line_fields_of_rna_conf_file

(* Parsing reference conf file *)
type ref_seq = {
  id : string ;
  species : string ;
  family : string ;
  path_fasta : string ;
}

let parse_line_fields_of_ref_seqs_file = function
  | [ id ; species ; family ; path_fasta ] ->
     { id ; species ; family ; path_fasta }
  | _ -> failwith "Syntax error in configuration file"

let parse_ref_seqs_file path =
  In_channel.read_lines path
  |> List.tl_exn
  |> List.map ~f:(String.split ~on:'\t')
  |> List.map ~f:parse_line_fields_of_ref_seqs_file
  

(* Run Trinity on RNA samples *)
let targets samples =
  let open Bistro_app in
  List.map samples ~f:(fun s ->
    let i = Bistro.Workflow.input s.path_fastq in
    let assembly = Trinity.trinity i in
    [
      [ "output" ; s.path_assembly ] %> assembly  ;
    ]
  )



(* RUN *)

(* Parsing rna conf file *)
let parsed_rna_conf_file = parse_rna_conf_file Sys.argv.(1)
(* Check if all paths exist*)
(* TO DO *)

(* Parsing reference sequences file *)
let parsed_ref_seqs_file = parse_ref_seqs_file Sys.argv.(2)
(* Check if all paths exist*)
(* TO DO *)


(* Run Trinity on RNA samples *)
let _ =
  parsed_rna_conf_file
  |> List.filter ~f:(fun s -> s.run_trinity )
  |> targets
  |> List.concat
  |> Bistro_app.simple

(* Run Seq_Dispatcher.py on Trinity assemblies *)
(** For each used reference species build a fasta with all its sequences **)
(*** Get each used reference species ***)
let used_ref_species = parsed_rna_conf_file
  |> List.map ~f:(fun s -> s.ref_species )
  |> List.sort_uniq compare
  
(* TO DO unique of a list*)

let _ = List.iter ~f:(printf "%s\n" ) used_ref_species





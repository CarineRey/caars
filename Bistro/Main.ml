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
let targets_trinity samples =
  let open Bistro_app in
  List.map samples ~f:(fun s ->
    let i = Bistro.Workflow.input s.path_fastq in
    let assembly = Trinity.trinity i in
    [
      [ "output" ; s.path_assembly ] %> assembly  ;
    ]
  )


(* cat fasta files *)
let cat_fasta_files (fasta_list: fasta workflow list) : fasta workflow = 
       workflow [
          cmd "cat" ~stdout:dest [ list ~sep:" " dep fasta_list  ] ;
    ]

let targets_cat fasta_list =
  let open Bistro_app in
    let cat = cat_fasta_files fasta_list in
    [
      [ "tmp" ; "ref_sequences" ; ] %> cat  ;
    ]

(* RUN *)

let check_path_exists path =
  if not (Sys.file_exists path) then (
    let msg = "No such file or directory " ^ path in
    failwith msg
  )

(* Parsing rna conf file *)
let parsed_rna_conf_file = parse_rna_conf_file Sys.argv.(1)

(* Check if all paths exist*)
let () =
  List.iter parsed_rna_conf_file ~f:(fun sample ->
      check_path_exists sample.path_fastq
    )

(* Parsing reference sequences file *)
let parsed_ref_seqs_file = parse_ref_seqs_file Sys.argv.(2)
(* Check if all paths exist*) (* TO DO *)

(** Get each used reference species **)
let used_ref_species = parsed_rna_conf_file
  |> List.map ~f:(fun s -> s.ref_species )
  |> Caml.List.sort_uniq compare
  
(** Concatenate all fasta with ref sequences **)
let all_ref_seqs_fasta = parsed_ref_seqs_file
  |> List.map ~f:(fun s -> s.path_fasta )
  |> Caml.List.sort_uniq compare

let fasta_paths = Caml.String.concat " " all_ref_seqs_fasta

let _ = let open Bistro_app in
   all_ref_seqs_fasta
	|> List.map ~f:(fun f -> Bistro.Workflow.input f)
	|> targets_cat
	|> Bistro_app.simple



(* Run Trinity on RNA samples *)
let _ =
  parsed_rna_conf_file
  |> List.filter ~f:(fun s -> s.run_trinity )
  |> targets_trinity
  |> List.concat
  |> Bistro_app.simple

(* Run Seq_Dispatcher.py on Trinity assemblies *)
(** For each used reference species build a fasta with all its sequences **)



let _ = List.iter ~f:(printf "%s\n" ) used_ref_species
let _ = print_string fasta_paths
let _ = List.iter ~f:(printf "\n%s\n" ) all_ref_seqs_fasta
(*let _ = print_string all_ref_seqs_fasta*)



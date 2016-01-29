open Core_kernel.Std
open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std

let trinity (fastq: _ fastq workflow) : fasta workflow =
	workflow [
		mkdir_p dest;
		cmd "Trinity" [ string "-single" ; dep fastq ; string "--seqType fq --max_memory 10G --full_cleanup " ; 
		string "--output" ; seq ~sep:"/" [ dest ; string "trinity"] ] ;
	]
    / selector [ "trinity.Trinity.fasta" ]


(* Parsing conf file *)
type sample = {
  id : string ;
  path_fastq : string ;
  run_trinity : bool ;
  path_assembly : string ;
}

let parse_line_fields_of_conf_file = function
  | [ id ; path_fastq ; run_trinity ; path_assembly ] ->
     let run_trinity = match run_trinity with
       | "yes" | "Yes" -> true
       | "no" | "No" -> false
       | _ -> failwith "Syntax error in configuration file" 
     in
     { id ; path_fastq ; run_trinity ; path_assembly }
  | _ -> failwith "Syntax error in configuration file"

let parse_conf_file path =
  In_channel.read_lines path
  |> List.tl_exn
  |> List.map ~f:(String.split ~on:'\t')
  |> List.map ~f:parse_line_fields_of_conf_file

    
let targets samples =
  let open Bistro_app in
  List.map samples ~f:(fun s ->
    let i = Bistro.Workflow.input s.path_fastq in
    let assembly = trinity i in
    let alignment = Mafft.mafft assembly in
    [
      [ "output" ; s.path_assembly ] %> assembly  ;
      [ "output" ; s.path_assembly ^ ".mafft" ] %> alignment ;
    ]
  )

let _ =
  parse_conf_file Sys.argv.(1)
  |> targets
  |> List.concat
  |> Bistro_app.simple


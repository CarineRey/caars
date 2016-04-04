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
  
(*
(* to run Trinity on RNA samples *)
let targets_trinity samples =
  let open Bistro_app in
  List.map samples ~f:(fun s ->
    let fastq = Bistro.Workflow.input s.path_fastq in
    let assembly = Trinity.trinity ~full_cleanup:true ~memory fastq in
    [
      [ "tmp_amalgam"; "assembly" ; s.path_assembly ] %> assembly  ;
    ]
  )

*)

(* to run makeblastdb on RNA samples *)

let parse_seqids = "yes"
let dbtype = "nucl"

let targets_apytram samples =
  let open Bistro_app in
  List.map samples ~f:(fun s ->
    let fastq = Bistro.Workflow.input s.path_fastq in
    let fasta = Trinity.fastool fastq in
    let db_blast = BlastPlus.makeblastdb dbtype s.species fasta in
    [
      [ "tmp_amalgam" ; "fasta"; s.species ^ ".fa" ] %> fasta  ;
      [ "tmp_amalgam" ; "BlastDB"; s.species ^ "_DB" ] %> db_blast  ;
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

  
(* Run ParseInput.py *)

module Amalgam = struct

type configuration = [ `configuration ] directory

let parse_input rna_conf_file species_tree_file ali_dir seq2sp_dir : configuration workflow= 
       workflow [
            cmd "../bin/ParseInput.py"  [ string rna_conf_file ;
                                          string species_tree_file;
                                          string ali_dir;
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

let seq_2_species_links =
    selector ["Validated_Sequences2Species"]

  
let norm_fasta_of_rna_conf rna_conf_file =
  parse_rna_conf_file rna_conf_file
  |> List.filter ~f:(fun s -> s.run_apytram || s.run_trinity)
  |> List.map ~f:(fun s ->
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
    

let seq_dispatcher ?tab_by_family_dir query query_species ref_transcriptome seq2fam : fasta workflow = 
       workflow [
            cmd "../bin/SeqDispatcher.py"  [ 
              option (opt "-tab_out_by_family" string) tab_by_family_dir;
			  opt "-q" string query ;
			  opt "-qs" string query_species ;
			  opt "-t" string ref_transcriptome ;
			  opt "-t2f" string seq2fam;
			  opt "-out" seq [ dest ; string ("/Trinity." ^ query_species )] ;
            ]
    ]

let trinity_annotated_fams_of_trinity_assemblies  =
    List.map ~f:(fun (s,trinity_assembly) ->
      let r = 
        seq_dispatcher 
          ~tab_by_family_dir:(configuration / seq_2_species_links) 
          trinity_assembly 
          s.species 
          target
          (configuration / 
      in
      (s, r)
  )
         

let main rna_conf_file species_tree_file ali_dir seq2sp_dir  =
    let configuration = parse_input rna_conf_file species_tree_file ali_dir seq2sp_dir in
    let norm_fasta = norm_fasta_of_rna_conf rna_conf_file in
    let trinity_assemblies = trinity_assemblies_of_norm_fasta norm_fasta in
    let open Bistro_app in 
    List.concat [
      [ [ "tmp_amalgam" ] %>  configuration ] ;
      List.map norm_fasta ~f:( fun (s,norm_fasta) ->
        [ "tmp_amalgam" ; "norm_fasta" ; s.species ^ ".fa" ] %> norm_fasta
        );
      List.map trinity_assemblies ~f:( fun (s,trinity_assemblies) ->
        [ "tmp_amalgam" ; "trinity_assemblies" ; "Trinity_assemblies." ^ s.species ^ ".fa" ] %> trinity_assemblies
        );
      
    
    ]
 

end

let  target_amalgam = Amalgam.main rna_conf_file species_tree_file ali_dir seq2sp_dir 

let _ = Bistro_app.local target_amalgam

(* Read Normalization of RNA sample

let _ = parsed_rna_conf_file
  |> List.filter ~f:(fun s -> s.run_apytram || s.run_trinity)
  |> targets_normalization
  |> List.concat
  |> Bistro_app.local
                                
*)

(* Run makeblastdb on RNA samples 
let _ =
  parsed_rna_conf_file
  |> List.filter ~f:(fun s -> s.run_apytram )
  |> targets_apytram
  |> List.concat
  |> Bistro_app.local


*)




(* Run Seq_Dispatcher.py on Trinity assemblies *)
(** For each used reference species build a fasta with all its sequences **)


(*let _ = List.iter ~f:(printf "%s\n" ) used_ref_species*)


(*let _ = print_string fasta_paths *)
(*let _ = List.iter ~f:(printf "\n%s\n" ) all_ref_seqs_fasta*)
(*let _ = print_string all_ref_seqs_fasta*)



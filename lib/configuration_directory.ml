open Core_kernel
open Bistro
open Bistro.Shell_dsl
open Wutils

type t = [`caars_configuration] directory

let make ~memory (config : Dataset.t) : t =
  let families_out = dest // "DetectedFamilies.txt" in
  let script = Bistro.Template_dsl.(
      [
        [seq ~sep:"\t" [string "Detected_families"; string "Fam_ID"]];
        List.map config.all_families ~f:(fun fam -> seq ~sep:"\t" [string fam.name; int fam.id ])
      ]
      |> List.concat
      |> seq ~sep:"\n"
    )
  in
  let ali_cmd_list = List.map config.all_families ~f:(fun fam ->
      let ali = Workflow.input (config.alignments_dir ^ "/" ^ fam.name ^ ".fa") in
      cmd "echo" [dep ali]
    )
  in
  let sp2seq_files = Sys.readdir config.seq2sp_dir
                     |> Array.filter ~f:(fun f ->
                         if Filename.check_suffix f ".tsv" then
                           true
                         else
                           false)
                     |> Array.to_list
  in
  let sp2seq_cmd_list = List.map sp2seq_files ~f:(fun sp2seq_file ->
      let sp2seq = Workflow.input (config.seq2sp_dir ^ "/" ^ sp2seq_file) in
      cmd "echo" [dep sp2seq]
    )
  in
  Workflow.shell ~np:1 ~descr:"parse_input" ~version:17 ~mem:(Workflow.int (memory * 1024)) (List.concat [
      [mkdir_p dest;
       cmd "python" ~img:caars_img [
         file_dump (string Scripts.parse_input) ;
         dep (Workflow.input config.sample_sheet) ;
         dep (Workflow.input config.species_tree_file) ;
         dep (Workflow.input config.alignments_dir) ;
         dep (Workflow.input config.seq2sp_dir) ;
         ident dest ;
       ];
       cmd "cp" [ file_dump script; families_out];
      ];
      ali_cmd_list;
      sp2seq_cmd_list;
    ])

let ref_transcriptome dir species =
  Workflow.select dir ["R_Sp_transcriptomes" ;  species ^ "_transcriptome.fa" ]

let ref_seq_fam_links dir species =
  Workflow.select dir ["R_Sp_Seq_Fam_links";  species ^ "_Fam_Seq.tsv"  ]

let ref_fams dir species family =
  Workflow.select dir ["R_Sp_Gene_Families"; species ^ "." ^ family ^ ".fa"]

let ali_species2seq_links dir family =
  Workflow.select dir ["Alignments_Species2Sequences" ; "alignments." ^  family ^ ".sp2seq.txt" ]

let family_metadata dir = Workflow.select dir ["FamilyMetadata.txt"]
let species_metadata dir = Workflow.select dir ["SpeciesMetadata.txt"]
let usable_families dir = Workflow.select dir ["UsableFamilies.txt"]
let detected_families dir = Workflow.select dir ["DetectedFamilies.txt"]

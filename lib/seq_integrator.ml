open Core_kernel
open Bistro
open Bistro.Shell_dsl
open Defs
open Wutils

let transform_species_list l =
  list ~sep:"," string l

let seq_integrator
    ?realign_ali
    ?resolve_polytomy
    ?species_to_refine_list
    ?no_merge
    ?merge_criterion
    ~family
    ~trinity_fam_results_dirs
    ~apytram_results_dir
    ~alignment_sp2seq
    alignment
  : [`seq_integrator] directory =

  let open Bistro.Shell_dsl in
  let merge_criterion_string =
    Option.bind merge_criterion ~f:(function
        | Merge -> None
        | Length -> Some "length"
        | Length_complete -> Some "length.complete"
      )
  in
  let get_trinity_file_list extension dirs =
    List.concat_map dirs ~f:(fun ((s : Rna_sample.t), dir) ->
        [ dep dir ; string ("/Trinity." ^ s.id ^ "." ^ s.species ^ "." ^ family ^ "." ^ extension) ; string ","]
      )
  in
  let get_apytram_file_list extension dir =
    [ dep dir ; string ("/apytram." ^ family ^ "." ^ extension) ; string ","]
  in
  let trinity_fasta_list = get_trinity_file_list "fa" trinity_fam_results_dirs in
  let trinity_sp2seq_list = get_trinity_file_list "sp2seq.txt" trinity_fam_results_dirs in

  let apytram_fasta = get_apytram_file_list "fa" apytram_results_dir in
  let apytram_sp2seq = get_apytram_file_list "sp2seq.txt" apytram_results_dir in

  let sp2seq = List.concat [[dep alignment_sp2seq ; string "," ] ; trinity_sp2seq_list ; apytram_sp2seq ]  in
  let fasta = List.concat [trinity_fasta_list; apytram_fasta]  in

  let tmp_merge = tmp // "tmp" in

  Workflow.shell ~version:12 ~descr:("SeqIntegrator.py:" ^ family) [
    mkdir_p tmp_merge ;
    cmd "python" ~img:caars_img [
      file_dump (string Scripts.seq_integrator);
      opt "-tmp" ident tmp_merge;
      opt "-log" seq [ tmp_merge ; string ("/SeqIntegrator." ^ family ^ ".log" )] ;
      opt "-ali" dep alignment ;
      opt "-fa" (seq ~sep:"") fasta;
      option (flag string "--realign_ali") realign_ali;
      option (opt "--merge_criterion" string) merge_criterion_string;
      option (flag string "--no_merge") no_merge;
      option (flag string "--resolve_polytomy") resolve_polytomy;
      opt "-sp2seq" (seq ~sep:"") sp2seq  ; (* list de sp2seq delimited by comas *)
      opt "-out" seq [ dest ; string "/" ; string family] ;
      option (opt "-sptorefine" transform_species_list) species_to_refine_list;
    ]
  ]

let seq_filter
    ?realign_ali
    ?resolve_polytomy
    ?species_to_refine_list
    ~filter_threshold
    ~family
    ~alignment
    ~tree
    ~sp2seq
  : [`seq_integrator] directory  =

  let open Bistro.Shell_dsl in
  let tmp_merge = tmp // "tmp" in

  Workflow.shell ~version:8 ~descr:("SeqFilter.py:" ^ family) [
    mkdir_p tmp_merge ;
    cmd "python" ~img:caars_img [
      file_dump (string Scripts.seq_filter);
      opt "-tmp" ident tmp_merge;
      opt "-log" seq [ tmp_merge ; string ("/SeqFilter." ^ family ^ ".log" )] ;
      opt "-ali" dep alignment ;
      opt "-t" dep tree;
      opt "--filter_threshold" float filter_threshold;
      option (flag string "--realign_ali") realign_ali;
      option (flag string "--resolve_polytomy") resolve_polytomy;
      opt "-sp2seq" dep sp2seq  ;
      opt "-out" seq [ dest ; string "/" ; string family] ;
      option (opt "-sptorefine" transform_species_list) species_to_refine_list;
    ]
  ]

let alignment dir (fam : Family.t) = Workflow.select dir [ fam.name ^ ".fa" ]
let tree dir (fam : Family.t) = Workflow.select dir [ fam.name ^ ".tree" ]
let sp2seq dir (fam : Family.t) = Workflow.select dir [ fam.name ^ ".sp2seq.txt" ]

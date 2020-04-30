open Core_kernel
open Bistro
open Wutils
    
let seq_dispatcher
    ?s2s_tab_by_family
    ~ref_db
    ~query
    ~query_species
    ~query_id
    ~ref_transcriptome
    ~threads
    ~seq2fam : [`seq_dispatcher] directory =
  let open Shell_dsl in
  Workflow.shell ~np:threads ~version:9 ~descr:("SeqDispatcher.py:" ^ query_id ^ "_" ^ query_species) [
    mkdir_p tmp;
    cmd "python" ~img:caars_img [
      file_dump (string Scripts.seq_dispatcher);
      option (flag string "--sp2seq_tab_out_by_family" ) s2s_tab_by_family;
      opt "-d" ident (seq ~sep:"," (List.map ref_db ~f:(fun blast_db -> seq [dep blast_db ; string "/db"]) ));
      opt "-tmp" ident tmp ;
      opt "-log" seq [ dest ; string ("/SeqDispatcher." ^ query_id ^ "." ^ query_species ^ ".log" )] ;
      opt "-q" dep query ;
      opt "-qs" string query_species ;
      opt "-qid" string query_id ;
      opt "-threads" ident np ;
      opt "-t" dep ref_transcriptome ;
      opt "-t2f" dep seq2fam;
      opt "-out" seq [ dest ; string ("/Trinity." ^ query_id ^ "." ^ query_species )] ;
    ]
  ]

let fasta_file_name (sample : Rna_sample.t) family =
  sprintf "Trinity.%s.%s.%s.fa" sample.id sample.species family

let get_fasta dir (sample : Rna_sample.t) family =
  Workflow.select dir [fasta_file_name sample family]

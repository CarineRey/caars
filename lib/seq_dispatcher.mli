open Bistro

val seq_dispatcher :
  ?s2s_tab_by_family:bool ->
  ref_db:'a Bistro.path Bistro.workflow list ->
  query:fasta file ->
  query_species:string ->
  query_id:string ->
  ref_transcriptome:fasta file ->
  threads:int ->
  seq2fam:'d Bistro.path Bistro.workflow ->
  [ `seq_dispatcher ] Bistro.directory

val fasta_file_name : Rna_sample.t -> string -> string
val get_fasta :  [ `seq_dispatcher ] Bistro.directory -> Rna_sample.t -> string -> fasta file

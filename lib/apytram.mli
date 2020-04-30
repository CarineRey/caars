open Bistro
open Wutils

type compressed_read_db = {
  s : Rna_sample.t ;
  concat_fasta : fasta file;
  index_concat_fasta : biopython_sequence_index file;
  rep_cluster_fasta : fasta file;
  reformated_cluster : fasta file;
  index_cluster : biopython_sequence_index file;
  cluster_rep_blast_db : blast_db file;
}

val apytram_multi_species :
  ?descr:string ->
  ?i:int ->
  ?evalue:float ->
  ?no_best_file:bool ->
  ?only_best_file:bool ->
  ?out_by_species:bool ->
  ?write_even_empty:bool ->
  ?id:float ->
  ?fid:float ->
  ?mal:int ->
  ?fmal:float ->
  ?len:float ->
  ?flen:float ->
  ?required_coverage:float ->
  ?stats:bool ->
  ?threads:int ->
  ?memory:int ->
  ?time_max:int ->
  query:fasta file ->
  fam:string ->
  compressed_read_db list ->
  [`apytram] directory

val get_fasta :
  [`apytram] directory ->
  family_name:string ->
  sample_id:string ->
  fasta file

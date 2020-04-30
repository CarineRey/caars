open Bistro

val fastq2fasta :
  ?descr:string ->
  ?dep_input:_ file -> (* FIXME: what's this for? *)
  fastq file ->
  fasta file

val fasta_read_normalization :
  ?descr:string ->
  max_cov:int ->
  threads:int ->
  ?memory:int ->
  ?max_memory:int ->
  fasta file OSE_or_PE.t ->
  fasta file OSE_or_PE.t

val trinity_fasta :
  ?tag:string ->
  ?full_cleanup:bool ->
  ?no_normalization:bool ->
  threads:int ->
  ?memory:int ->
  fasta file OSE_or_PE.t ->
  fasta file

val assembly_stats :
  ?descr:string ->
  fasta file ->
  text file

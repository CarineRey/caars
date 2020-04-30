(** Utility functions to build workflows *)

open Bistro

class type blast_db = object
  inherit text
  method format : [`blast_db]
end

class type biopython_sequence_index = object
  inherit binary_file
  method format : [`biopython_sequence_index]
end

val caars_img : Shell_dsl.container_image list
val descr : ?tag:string -> string -> string

val fasta_concat : ?tag:string -> fasta file list -> fasta file
val fastq_concat : ?tag:string -> fastq file list -> fastq file

val build_biopythonindex :
  ?tag:string ->
  fasta file ->
  biopython_sequence_index file

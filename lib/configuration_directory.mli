(**
   Early step in analysis, to generate many files that subsequent steps take as input
*)

open Bistro

type t = [`caars_configuration] directory

val make : memory:int -> Dataset.t -> t

val ref_transcriptome : t -> string -> fasta file
val ref_seq_fam_links : t -> string -> fasta file
val ref_fams : t -> string -> string -> fasta file
val ali_species2seq_links : t -> string -> fasta file

val usable_families : t -> text file
val family_metadata : t -> text file
val species_metadata : t -> text file
val detected_families : t -> text file


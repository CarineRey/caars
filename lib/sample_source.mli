(** Description of a data source for sequencing data *)

open Defs

type t =
  | Fastq_file of string OSE_or_PE.t
  | Fasta_file of string OSE_or_PE.t

val orientation : t -> (single_end_orientation, paired_end_orientation) either

type t =
  | Fastq_file of string OSE_or_PE.t
  | Fasta_file of string OSE_or_PE.t

let orientation = function
  | Fastq_file x
  | Fasta_file x -> OSE_or_PE.orientation x

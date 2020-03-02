open Bistro
open Commons

val alignement_fasta : string -> output directory -> fasta file
val gene_tree : string -> output directory -> newick file
val sp2seq_link : string -> output directory -> sp2seq_link file
val build_term : Configuration.t -> unit workflow

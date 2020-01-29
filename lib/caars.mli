open Bistro
open Commons

val alignement_fasta : string -> output dworkflow -> fasta pworkflow
val gene_tree : string -> output dworkflow -> [`newick] pworkflow
val sp2seq_link : string -> output dworkflow -> sp2seq_link pworkflow
val build_term : Configuration.t -> unit workflow

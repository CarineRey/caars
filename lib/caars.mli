open Bistro.Std
open Bistro_bioinfo.Std
open Bistro_utils
open Commons

val alignement_fasta : string -> (output, fasta) selector
val gene_tree : string -> (output, [`newick]) selector
val sp2seq_link : string -> (output, sp2seq_link) selector
val build_app : Configuration.t -> unit Bistro_app.t

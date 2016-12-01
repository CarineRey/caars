open Bistro.Std
open Bistro_bioinfo.Std

type output = [ `amalgam_output ]

val alignement_fasta : string -> (output, fasta) selector

val gene_tree : string -> (output, [`newick]) selector

type sp2seq_link
val sp2seq_link : string -> (output, sp2seq_link) selector


val build_app : Configuration.t -> unit Bistro_app.t

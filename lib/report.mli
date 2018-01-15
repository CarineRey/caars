open Commons
open Bistro_utils

val generate :
  trinity_assemblies_stats:(rna_sample * Trinity.assembly_stats Term.path) list ->
  final_plots:_ Term.path ->
  string ->
  unit


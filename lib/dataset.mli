open Defs

type t = {
  samples : Rna_sample.t list ;
  sample_sheet : string ;
  alignments_dir : string ;
  seq2sp_dir : string ;
  species_tree_file : string ;
  reference_species : string list ;
  all_families : Family.t list ;
  used_families : Family.t list ;
}

val apytram_samples : t -> Rna_sample.t list
val apytram_groups : t -> string list
val apytram_reference_species : t -> string list list

val trinity_samples : t -> Rna_sample.t list
val reference_samples : t -> Rna_sample.t list
val has_at_least_one_sample_with_reference : t -> bool

val make :
  ?family_subset_file:string ->
  sample_sheet:string ->
  species_tree_file:string ->
  alignments_dir:string ->
  seq2sp_dir:string ->
  unit -> t

open Defs
open Bistro

val seq_integrator :
  ?realign_ali:bool ->
  ?resolve_polytomy:bool ->
  ?species_to_refine_list:string list ->
  ?no_merge:bool ->
  ?merge_criterion:merge_criterion ->
  family:string ->
  trinity_fam_results_dirs:(Rna_sample.t * [`seq_dispatcher] directory) list ->
  apytram_results_dir:fasta file ->
  alignment_sp2seq:fasta file ->
  fasta file ->
  [ `seq_integrator ] directory

val seq_filter :
  ?realign_ali:bool ->
  ?resolve_polytomy:bool ->
  ?species_to_refine_list:string list ->
  filter_threshold:float ->
  family:string ->
  alignment:fasta file ->
  tree:text file ->
  sp2seq:text file ->
  [ `seq_integrator ] directory

val alignment : [ `seq_integrator ] directory -> Family.t -> fasta file
val tree : [ `seq_integrator ] directory -> Family.t -> text file
val sp2seq : [ `seq_integrator ] directory -> Family.t -> text file

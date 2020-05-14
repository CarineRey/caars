open Bistro

val just_parse_workflow :
  outdir:string ->
  Pipeline.t ->
  unit workflow

val full_analysis_workflow :
  ?get_reads:bool ->
  ?debug:bool ->
  outdir:string ->
  Pipeline.t ->
  unit workflow

val main :
  sample_sheet:string ->
  outdir:string ->
  species_tree_file:string ->
  alignments_dir:string ->
  seq2sp_dir:string ->
  ?np:int ->
  ?memory:int ->
  no_reconcile:bool ->
  ?ali_sister_threshold:float ->
  ?merge_criterion:string ->
  debug:bool ->
  get_reads:bool ->
  just_parse_input:bool ->
  ?html_report:string ->
  quiet:bool ->
  use_docker:bool ->
  ?family_to_use:string ->
  unit -> unit

val command : Core.Command.t

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

val command : Core.Command.t

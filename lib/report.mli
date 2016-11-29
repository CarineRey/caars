open Configuration

val generate :
  trinity_assemblies_stats:(rna_sample * Trinity.assembly_stats Bistro_app.path) list ->
  string ->
  unit


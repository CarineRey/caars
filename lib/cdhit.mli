open Bistro

val lap :
  ?tag:string ->
  ?p:float ->
  ?m:float ->
  ?d:float ->
  fasta file ->
  [`cdhit_lap] directory

val cluster_rep_of_lap : [`cdhit_lap] directory -> fasta file
val cluster_of_lap : [`cdhit_lap] directory -> text file

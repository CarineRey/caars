type t = {
  n50 : int option ;
  nb_genes : int option ;
  nb_transcripts : int option ;
  gc : float option ;
}
val parse : string -> t

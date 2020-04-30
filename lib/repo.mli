val just_parse_repo :
  Pipeline.t ->
  Bistro_utils.Repo.t

val full_repo :
  ?get_reads:bool ->
  ?debug:bool ->
  Pipeline.t -> 
  Bistro_utils.Repo.t

val precious_repo :
  Pipeline.t -> 
  Bistro_utils.Repo.t
  

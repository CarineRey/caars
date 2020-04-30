type entry = {
  id : string ;
  descr : string ;
  start_time : float ;
  end_time : float ;
  size : int ;
}
[@@deriving sexp]

val create : string -> Bistro_engine.Logger.t

val read_db : string -> entry list

open Base

type entry = {
  id : string ;
  descr : string ;
  start_time : float ;
  end_time : float ;
  size : int ;
}
[@@deriving sexp]

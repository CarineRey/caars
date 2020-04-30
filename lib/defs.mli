(** Various basic definitions *)

type url = string

type ('a,'b) either =
  | Left of 'a
  | Right of 'b

type single_end_orientation =
  | F
  | R
  | US

type paired_end_orientation =
  | FR
  | RF
  | UP

module Family : sig
  type t = {
    name : string ;
    id : int;
  }
end

type merge_criterion =
  | Merge
  | Length
  | Length_complete

val merge_criterion_of_string : string -> merge_criterion option

(** concatenates strings using underscore as separator *)
val id_concat : string list -> string

(** {5 Association lists} *)

type 'a assoc = (string * 'a) list

val assoc : 'a list -> f:('a -> 'b) -> ('a * 'b) list
val assoc_opt : 'a list -> f:('a -> 'b option) -> ('a * 'b) list
val assoc_map : ('a * 'b) list -> f:('a -> 'b -> 'c) -> ('a * 'c) list
val assoc_filter_map : ('a * 'b) list -> f:('a -> 'b -> 'c option) -> ('a * 'c) list
val ( $ ) : ('a * 'b) list -> 'a -> 'b

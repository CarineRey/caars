(** Values that are oriented and that can be single or paired end *)

open Defs

type 'a t =
  | Single_end of {
      reads : 'a ;
      orientation : single_end_orientation ;
    }
  | Paired_end of {
      reads1 : 'a ;
      reads2 : 'a ;
      orientation : paired_end_orientation ;
    }

val se : 'a -> single_end_orientation -> 'a t
val pe : 'a -> 'a -> paired_end_orientation -> 'a t

val map : 'a t -> f:('a -> 'b) -> 'b t
val orientation : _ t -> (single_end_orientation, paired_end_orientation) either

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

let se reads orientation = Single_end { reads ; orientation }
let pe reads1 reads2 orientation = Paired_end { reads1 ; reads2 ; orientation }
let map x ~f = match x with
  | Single_end se -> Single_end { se with reads = f se.reads }
  | Paired_end pe -> Paired_end { pe with reads1 = f pe.reads1 ; reads2 = f pe.reads2 }

let orientation = function
  | Single_end se -> Left se.orientation
  | Paired_end pe -> Right pe.orientation

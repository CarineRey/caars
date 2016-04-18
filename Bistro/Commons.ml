type ('a,'b) either =
  | Left of 'a
  | Right of 'b

type orientation_single =
  | F
  | R
  | US

type orientation_paired =
  | FR
  | RF
  | UP

type 'a sample_fastq =
  | Single_end of 'a * orientation_single
  | Paired_end of 'a * 'a * orientation_paired

let sample_fastq_map f = function 
  | Single_end ( x , o ) ->  Single_end ( f x , o )
  | Paired_end ( lx, rx, o ) -> Paired_end (f lx, f rx,	o)

let sample_fastq_is_paired = function
  | Single_end _ ->  false
  | Paired_end _ -> true

let sample_fastq_orientation = function
  | Single_end ( x , o ) ->  Left o 
  | Paired_end ( lx, rx, o ) -> Right o

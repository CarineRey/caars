open Core_kernel

type url = string

type 'a assoc = (string * 'a) list

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

module Family = struct
  type t = {
    name : string ;
    id : int;
  }
end

type merge_criterion =
  | Merge
  | Length
  | Length_complete

let merge_criterion_of_string = function
  | "merge" -> Some Merge
  | "length" -> Some Length
  | "length.complete" -> Some Length_complete
  | _ -> None

let ( $ ) a k = List.Assoc.find_exn ~equal:Poly.equal a k

let assoc keys ~f =
  List.map keys ~f:(fun k -> k, f k)

let assoc_map a ~f =
  List.map a ~f:(fun (k, v) -> k, f k v)

let assoc_filter_map a ~f =
  List.filter_map a ~f:(fun (k, v) -> Option.map (f k v) ~f:(fun v -> k, v))

let assoc_opt keys ~f =
  List.filter_map keys ~f:(fun k ->
      Option.map (f k) ~f:(fun v -> k, v)
    )

let id_concat xs =
  List.filter xs ~f:(String.( <> ) "")
  |> String.concat ~sep:"_"


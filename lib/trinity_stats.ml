open Core

type t = {
  n50 : int option ;
  nb_genes : int option ;
  nb_transcripts : int option ;
  gc : float option ;
}

let items = [
  `Total_trinity_genes, "Total trinity 'genes'", `int;
  `N50, "Contig N50", `int;
  `Total_trinity_transcripts, "Total trinity transcripts", `int;
  `GC, "Percent GC", `float;
 ]

let parse_int s = `int (Int.of_string s)
let parse_float s = `float (Float.of_string s)
let parse_value value_type s =
  let f = match value_type with
    | `int -> parse_int
    | `float -> parse_float
  in
  f (String.strip s)

let extract_item lines (key, marker, value_type) =
  List.find_map lines ~f:(fun l ->
      match String.lsplit2 l ~on:':' with
      | None -> None
      | Some (left, right) ->
        if String.(String.strip left = marker)
        then Some (key, parse_value value_type right)
        else None
    )


let parse_aux items lines =
  List.filter_map items ~f:(extract_item lines)

let find_int dict key =
  match List.Assoc.find ~equal:Poly.( = ) dict key with
  | None -> None
  | Some (`float _) -> None
  | Some (`int i) -> Some i

let find_float dict key =
  match List.Assoc.find ~equal:Poly.( = ) dict key with
  | None -> None
  | Some (`float f) -> Some f
  | Some (`int _) -> None


let parse fn =
  let lines = In_channel.read_lines fn in
  let dict = parse_aux items lines in
  {
    n50 = find_int dict `N50 ;
    nb_genes = find_int dict `Total_trinity_genes ;
    nb_transcripts = find_int dict `Total_trinity_transcripts ;
    gc = find_float dict `GC ;
  }

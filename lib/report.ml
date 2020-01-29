open Core
open Tyxml_html

let k = txt

let optint = function
  | None -> k"NA"
  | Some i -> k (Int.to_string i)

let optfloat = function
  | None -> k"NA"
  | Some i -> k (Float.to_string i)

(* let svg_from_file fn =
 *   let contents = In_channel.read_all fn in
 *   let src = sprintf "data:image/%s;base64,%s" "svg+xml" contents in
 *   Tyxml_html.img ~src ~alt:"" () *)

let head t =
  head (title (txt t)) [
    link ~rel:[`Stylesheet] ~href:"http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap.min.css" () ;
    link ~rel:[`Stylesheet] ~href:"http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap-theme.min.css" () ;
    script ~a:[a_src "https://code.jquery.com/jquery.js"] (txt "") ;
    script ~a:[a_src "http://netdna.bootstrapcdn.com/bootstrap/3.0.2/js/bootstrap.min.js"] (txt "") ;
  ]

let trinity_section trinity_assemblies_stats =
  let table_headers = [
    tr [
      th [ k "Species" ] ;
      th [ k "nb_genes" ] ;
      th [ k "nb_transcripts" ] ;
      th [ k "gc"] ;
      th [ k "n50" ] ;
    ]
  ]
  in
  let foreach_sample (sample, assembly_stats) =
    let { Trinity_stats.n50 ; nb_genes; gc; nb_transcripts } = Trinity_stats.parse assembly_stats in
    tr [
      td [ k sample.Commons.species ] ;
      td [ optint nb_genes ] ;
      td [ optint nb_transcripts ] ;
      td [ optfloat gc] ;
      td [ optint n50 ] ;
    ]
  in
  [
    h2 [ k"Trinity assemblies stats" ] ;
    table ~a:[a_class ["table" ; "table-condensed"]] (List.concat [ table_headers; (List.map trinity_assemblies_stats ~f:foreach_sample)])
  ]

(* http://ocsigen.org/tyxml/4.0.1/manual/intro*)
let render ~trinity_assemblies_stats =
  let mytitle = "Caars report" in
  let contents = List.concat [
      [
        h1 [txt "A fabulous title"] ;
        txt "This is a fabulous content." ;
      ] ;
      trinity_section trinity_assemblies_stats ;
      
    ]

    in
  html
    (head mytitle)
    (body [ div ~a:[a_class ["container"]] contents ])



let save path doc =
  let buf = Buffer.create 253 in
  let formatter = Format.formatter_of_buffer buf in
  Tyxml_html.pp () formatter doc ;
  Out_channel.with_file path ~f:(fun oc ->
      let contents = Buffer.contents buf in
      Out_channel.output_string oc contents
    )

let generate ~trinity_assemblies_stats dest =
  let doc = render ~trinity_assemblies_stats in
  save dest doc

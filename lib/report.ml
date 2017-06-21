open Core.Std
open Tyxml_html
open Commons

let k = pcdata

let head t =
  head (title (pcdata t)) [
    link ~rel:[`Stylesheet] ~href:"http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap.min.css" () ;
    link ~rel:[`Stylesheet] ~href:"http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap-theme.min.css" () ;
    script ~a:[a_src "https://code.jquery.com/jquery.js"] (pcdata "") ;
    script ~a:[a_src "http://netdna.bootstrapcdn.com/bootstrap/3.0.2/js/bootstrap.min.js"] (pcdata "") ;
  ]

let trinity_section trinity_assemblies_stats =
  let optint = function
    | None -> k"NA"
    | Some i -> k (Int.to_string i)
  in
  let optfloat = function
    | None -> k"NA"
    | Some i -> k (Float.to_string i)
  in
  let table_headers =
  [tr [
      th [ k "Species" ] ;
      th [ k "nb_genes" ] ;
      th [ k "nb_transcripts" ] ;
      th [ k "gc"] ;
      th [ k "n50" ] ;
    ]
    ]
  in
  let foreach_sample (sample, Bistro_app.Path assembly_stats) =
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
  let contents =
    List.concat [
      [
        h1 [pcdata "A fabulous title"] ;
        pcdata "This is a fabulous content." ;
      ] ;
      trinity_section trinity_assemblies_stats ;
    ]
    |> div ~a:[a_class ["container"]]

    in
  html
    (head mytitle)
    (body [ contents ])



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



open Core.Std
open Tyxml_html

let k = pcdata

let head t =
  head (title (pcdata t)) [
    link ~rel:[`Stylesheet] ~href:"http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap.min.css" () ;
    link ~rel:[`Stylesheet] ~href:"http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap-theme.min.css" () ;
    script ~a:[a_src "https://code.jquery.com/jquery.js"] (pcdata "") ;
    script ~a:[a_src "http://netdna.bootstrapcdn.com/bootstrap/3.0.2/js/bootstrap.min.js"] (pcdata "") ;
  ]

let trinity_section trinity_assemblies_stats =
  let foreach_sample (sample, trinity_output) =
    tr [
      td [ k sample.Configuration.species ] ;
      td [ k "42" ] ;
    ]
  in
  [
    h2 [ k"Trinity assemblies stats" ] ;
    table ~a:[a_class ["table"]] (List.map trinity_assemblies_stats ~f:foreach_sample)
  ]

(* http://ocsigen.org/tyxml/4.0.1/manual/intro*) 
let render ~trinity_assemblies_stats =
  let mytitle = "Amalgam report" in
  let contents =
    List.concat [
      [
        h1 [pcdata "A fabulous title"] ;
        pcdata "This is a fabulous content." ;
      ] ;
      trinity_section trinity_assemblies_stats ;
    ]
    |> div ~a:[a_class ["content"]]

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



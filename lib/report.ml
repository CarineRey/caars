open Core.Std
open Tyxml_html

let head t =
  head (title (pcdata t)) [
    link ~rel:[`Stylesheet] ~href:"http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap.min.css" () ;
    link ~rel:[`Stylesheet] ~href:"http://netdna.bootstrapcdn.com/bootstrap/3.0.2/css/bootstrap-theme.min.css" () ;
    script ~a:[a_src "https://code.jquery.com/jquery.js"] (pcdata "") ;
    script ~a:[a_src "http://netdna.bootstrapcdn.com/bootstrap/3.0.2/js/bootstrap.min.js"] (pcdata "") ;
  ]

let render ~stats =
  let contents = [] in
  html
    (head "Amalgam report")
    (body [ div ~a:[a_class ["container"]] contents ])


let save path doc =
  let buf = Buffer.create 253 in
  let formatter = Format.formatter_of_buffer buf in
  Tyxml_html.pp () formatter doc ;
  Out_channel.with_file path ~f:(fun oc ->
      let contents = Buffer.contents buf in
      Out_channel.output_string oc contents
    )

let generate ~stats dest =
  let doc = render ~stats in
  save dest doc



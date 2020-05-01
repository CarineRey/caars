open Core_kernel
open Tyxml_html

let save doc dest =
  let buf = Buffer.create 253 in
  let formatter = Format.formatter_of_buffer buf in
  Tyxml_html.pp () formatter doc ;
  Out_channel.write_all dest ~data:(Buffer.contents buf)

let encode_trace trace =
  [%sexp_of: Caars_execution_trace.entry list] trace
  |> Sexp.to_string_mach

let generate trace dest =
  let head =
    head (title (txt "Caars execution report")) [
      script (cdata_script Js_code.string) ;
      link ~rel:[`Stylesheet] ~href:"https://stackpath.bootstrapcdn.com/bootstrap/3.4.1/css/bootstrap.min.css" () ;
      script ~a:[a_src "https://code.jquery.com/jquery-1.12.4.min.js"] (txt "") ;
      script ~a:[a_src "https://stackpath.bootstrapcdn.com/bootstrap/4.4.1/js/bootstrap.min.js"] (txt "") ;
    ]
  in
  let data = encode_trace trace in
  let doc = html head (body ~a:[a_user_data "events" data] [div ~a:[a_id "container" ; a_class ["container"]] []]) in
  save doc dest

open Base
open Js_browser

type trace = Caars_execution_trace.entry list

let sp = Printf.sprintf
let float = Float.of_int

let elapsed_time (e : Caars_execution_trace.entry) =
  e.end_time -. e.start_time

let shorten_id id =
  if String.length id <=  7 then id
  else String.sub id ~pos:0 ~len:7

let div_str x y =
  Float.to_string_hum ~delimiter:' ' ~strip_zero:true ~decimals:1 (x /. y)

let format_duration d =
  let minute = 60. in
  let hour = 60. *. minute in
  let day = 24. *. hour in
  if Float.(d > day) then sp "%sd" (div_str d day)
  else if Float.(d > hour) then sp "%sh" (div_str d hour)
  else if Float.(d > minute) then sp "%smin" (div_str d minute)
  else sp "%.0fs" d

let format_size s =
  let ko = 1024 in
  let mo = ko * 1024 in
  let go = mo * 1024 in
  if s > go then sp "%sGB" (div_str (float s) (float go))
  else if s > mo then sp "%sMB" (div_str (float s) (float mo))
  else if s > ko then sp "%sKB" (div_str (float s) (float ko))
  else sp "%dB" s

module Task_table = struct
  open Vdom

  let cols = ["id" ; "descr" ; "elapsed time" ; "size"]

  let table_head =
    elt "thead" [
      elt "tr" (List.map cols ~f:(fun s -> elt "th" [text s]))
    ]

  let string_td s = elt "td" [ text s ]

  let table_line (e : Caars_execution_trace.entry) =
    elt "tr" [
      string_td (shorten_id e.id) ;
      string_td e.descr ;
      string_td (format_duration (elapsed_time e)) ;
      string_td (format_size e.size) ;
    ]
      
  let view data =
    let lines = List.map data ~f:table_line in
    elt "table" ~a:[attr "class" "table"] [
      table_head ;
      elt "tbody" lines
    ]
end

module Selection_stats = struct
  let view data =
    let n = List.length data in
    let total_time = List.fold data ~init:0. ~f:(fun acc e -> acc +. elapsed_time e) in
    let total_space = List.fold data ~init:0 ~f:(fun acc e -> acc + e.size) in
    let open Vdom in
    elt "table" ~a:[attr "class" "table .col-md-8 .table-bordered .table-condensed"] [
      elt "tbody" [
        elt "tr" [ elt "th" [text "Number of steps"] ;
                   elt "td" [text (Int.to_string n) ] ] ;
        elt "tr" [ elt "th" [text "Total time"] ;
                   elt "td" [text (format_duration total_time)] ] ;
        elt "tr" [ elt "th" [text "Total space"] ;
                   elt "td" [text (format_size total_space)] ] ;
      ]
    ]
end

type model = {
  data : trace ;
  selection : string ;
}

let filter_data (data : trace) sel =
  let res =
    String.split ~on:',' sel
    |> List.map ~f:(fun sel -> Re.Glob.glob sel |> Re.compile)
  in
  List.filter data ~f:(fun d ->
      List.exists res ~f:(fun re -> Re.execp re d.descr)
    )

let filter_form selection =
  let open Vdom in
  elt "form" ~a:[attr "class" "form-inline"] [
    div ~a:[attr "class" "form-group"] [
      elt "label" ~a:[attr "for" "selection-input"] [text "Filter"] ;
      input [] ~a:[attr "class" "form-control" ;
                   attr "id" "selection-input" ;
                   oninput (fun s -> `Selection_input s);
                   value selection ] ;
    ]
  ]

let view { data ; selection } =
  let open Vdom in
  let selected_data = filter_data data selection in
  div [
    elt "h1" [text "CAARS execution report"] ;
    div ~a:[class_ "row align-items-end"] [
      div ~a:[class_ "col-md-4"] [ filter_form selection ];
      div ~a:[class_ "col-md-8"] [ Selection_stats.view selected_data ];
    ] ;
    Task_table.view selected_data ;
  ]

let update m = function
  | `Selection_input selection -> { m with selection }

let run () =
  let data =
    Element.get_attribute (Document.body document) "data-events"
    |> Parsexp.Single.parse_string_exn
    |> [%of_sexp: Caars_execution_trace.entry list]
  in
  let init = { data ; selection = "" } in
  let app = Vdom.simple_app ~init ~view ~update () in
  match Document.get_element_by_id document "container" with
  | Some container ->
    Vdom_blit.run app
    |> Vdom_blit.dom
    |> Element.append_child container
  | None -> Window.alert window "container elt not found!"

let () = Window.set_onload window run

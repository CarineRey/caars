open Core_kernel
open Bistro_engine

type entry = {
  id : string ;
  descr : string ;
  start_time : float ;
  end_time : float ;
  size : int ;
}
[@@deriving sexp]

let save_entry fn e =
  let handler = Dbm.(opendbm fn [Dbm_rdwr ; Dbm_create]) 0o600 in
  Exn.protect ~f:(fun () ->
      let s = sexp_of_entry e |> Sexp.to_string_mach in
      Dbm.replace handler e.id s
    )
    ~finally:(fun () -> Dbm.close handler)

let sizeof db id =
  match Bistro_engine.Misc.du (Db.cache db id) with
  | Ok i -> i
  | Error (`Msg msg) -> failwith msg

let entry_of_task_result db start _end_ : Task_result.t -> entry option = function
  | Shell s ->
    {
      id = s.id ;
      descr = s.descr ;
      start_time = start ;
      end_time = _end_ ;
      size = sizeof db s.id ;
    }
    |> Option.some
  | Plugin s ->
    {
      id = s.id ;
      descr = s.descr ;
      start_time = start ;
      end_time = _end_ ;
      size = sizeof db s.id ;
    }
    |> Option.some

  | Input _
  | Select _
  | Container_image_fetch _ -> None

let create fn =
  object
  method event db _ = function
    | Logger.Workflow_ended { outcome ; start ; _end_ } ->
      let maybe_entry = entry_of_task_result db start _end_ outcome in
      Option.iter maybe_entry ~f:(save_entry fn)
    | Workflow_started _
    | Workflow_ready _
    | Workflow_skipped _
    | Workflow_allocation_error _
    | Workflow_collected _
    | Debug _
    | Singularity_image_collected _ -> ()

  method stop = Lwt.return ()
end

let read_db fn =
  let handler = Dbm.(opendbm fn [Dbm_rdwr]) 0o600 in
  Exn.protect ~f:(fun () ->
      let acc = ref [] in
      let f _ data =
        let e = Sexp.of_string data |> entry_of_sexp in
        acc := e :: !acc
      in
      Dbm.iter f handler ;
      !acc
    )
    ~finally:(fun () -> Dbm.close handler)

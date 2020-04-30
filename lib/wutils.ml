open Core_kernel
open Bistro

class type blast_db = object
  inherit text
  method format : [`blast_db]
end

let caars_img = Bistro.Shell_dsl.[ docker_image ~account:"carinerey" ~name:"caars_env" ~tag:"master_20200421" () ]

let descr ?tag d =
  match tag with
  | None -> d
  | Some tag -> sprintf "%s:%s" d tag

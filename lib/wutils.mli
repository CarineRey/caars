(** Utility functions to build workflows *)

open Bistro

class type blast_db = object
  inherit text
  method format : [`blast_db]
end

val caars_img : Shell_dsl.container_image list
val descr : ?tag:string -> string -> string

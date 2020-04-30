open Core_kernel
open Bistro
open Bistro.Shell_dsl

class type blast_db = object
  inherit text
  method format : [`blast_db]
end

class type biopython_sequence_index = object
  inherit binary_file
  method format : [`biopython_sequence_index]
end

let caars_img = Bistro.Shell_dsl.[ docker_image ~account:"carinerey" ~name:"caars_env" ~tag:"master_20200421" () ]

let descr ?tag d =
  match tag with
  | None -> d
  | Some tag -> sprintf "%s:%s" d tag

let concat ?tag xs =
  let descr = descr ?tag "concat" in
  match xs with
  | [] -> Workflow.shell ~descr [ cmd "touch" [ dest ] ]
  | x :: [] -> x
  | fXs ->
    Workflow.shell ~descr [
      cmd "cat" ~stdout:dest [ list dep ~sep:" " fXs ]
    ]

let fasta_concat = concat
let fastq_concat = concat

(* Fasta file and its index must be in the same directory due to
   biopython wich retains the relative path between these 2 files.
   A different location is incompatible with the bistro docker usage
   workflow by worflow. To avoid to cp the complete fasta file we
   use a symbolic link. *)
let build_biopythonindex ?tag (fasta : fasta file) : biopython_sequence_index file =
  Workflow.shell ~version:1 ~descr:(descr ?tag "build_biopythonindex_fasta.py") [
    mkdir_p dest ;
    within_container caars_img (
      and_list [
        cmd "ln" [ string "-s" ; dep fasta ; dest // "seq.fa" ] ;
        cmd "python" [
          file_dump (string Scripts.build_biopythonindex_fasta) ;
          dest // "index" ;
          dest // "seq.fa" ]
      ]
    )
  ]

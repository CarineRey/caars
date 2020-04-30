open Core_kernel
open Bistro
open Bistro.Shell_dsl
open Wutils

let lap ?tag ?p ?m ?d (fasta : fasta file) =
  let out = dest // "cluster_rep.fa" in
  Workflow.shell ~version:1 ~descr:(descr ?tag "cdhit-lap") [
    mkdir_p dest;
    cmd "cd-hit-lap" ~img:caars_img [
      opt "-i" dep fasta;
      opt "-o" Fn.id out ;
      option ( opt "-p" float ) p;
      option ( opt "-m" float ) m;
      option ( opt "-d" float ) d;
    ]
  ]

let cluster_rep_of_lap dir = Workflow.select dir ["cluster_rep.fa"]
let cluster_of_lap dir = Workflow.select dir ["cluster_rep.fa.clstr"]

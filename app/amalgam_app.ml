(*
# File: amalgam.ml
# Created by: Carine Rey
# Created on: March 2016
#
#
# Copyright 2016 Carine Rey
# This software is a computer program whose purpose is to assembly
# sequences from RNA-Seq data (paired-end or single-end) using one or
# more reference homologous sequences.
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
*)

open Core.Std
open Bistro.Std
open Bistro.EDSL
open Bistro_bioinfo.Std
open Commons

let main sample_sheet outdir species_tree_file alignments_dir seq2sp_dir np memory no_reconcile refinetree (*refineali*) ali_sister_threshold debug just_parse_input html_report dag_dot quiet () =
  let logger quiet html_report dag_dot =
  Bistro_logger.tee
    (if quiet then Bistro_logger.null else Bistro_console_logger.create ())
    (Bistro_logger.tee
       (match dag_dot with
        | Some path -> Bistro_dot_output.create path
        | None -> Bistro_logger.null)
       (match html_report with
        | Some path -> Bistro_html_logger.create path
        | None -> Bistro_logger.null)
    )
  in
  let run_reconciliation = match no_reconcile with
    | true -> false
    | false -> true
  in
  (*let refinetree = Option.value ~default:false refinetree in
  let refineali = Option.value ~default:false refineali in
  let debug = Option.value ~default:false debug in
  let just_parse_input = Option.value ~default:false just_parse_input in *)
  let refineali = false in
  let ali_sister_threshold = Option.value ~default:0.0 ali_sister_threshold in
  let np = Option.value ~default:2 np in
  let memory = Option.value ~default:1 memory in
  let configuration = Configuration.load ~sample_sheet ~species_tree_file ~alignments_dir ~seq2sp_dir ~np ~memory ~run_reconciliation ~debug ~just_parse_input ~refinetree ~refineali ~ali_sister_threshold ~outdir in
  let amalgam_app = Amalgam.build_app configuration in
  Bistro_app.(
    run ~logger:(logger quiet html_report dag_dot) ~np:configuration.Configuration.threads ~mem:(1024 * configuration.Configuration.memory) ~keep_all:false ~bistro_dir:"_AmalgamCache" amalgam_app
  )

let spec =
  let open Command.Spec in
  empty
  +> flag "--sample-sheet"    (required file)   ~doc:"PATH sample sheet file."
  +> flag "--outdir"          (required string) ~doc:"PATH Destination directory."
  +> flag "--species-tree"    (required file)   ~doc:"ABSOLUTE PATH Species tree file in nw format containing all species. Warning absolute path is required."
  +> flag "--alignment-dir"   (required string) ~doc:"PATH Directory containing all gene family alignments (Family_name.fa) in fasta format."
  +> flag "--seq2sp-dir"      (required string) ~doc:"PATH Directory containing all link files (Family_name.tsv). A line for each sequence and its species spaced by a tabulation."
  +> flag "--np"              (optional int)    ~doc:"INT Number of CPUs (at least 2). (Default:2)"
  +> flag "--memory"          (optional int)    ~doc:"INT Number of GB of system memory to use.(Default:1)"
  +> flag "--no-reconcile"    no_arg            ~doc:"Not run final Reconciliation step"
  +> flag "--refinetree"      no_arg            ~doc:"Refine topology during final Reconciliation step (Default:false)"
(*  +> flag "--refineali"       no_arg            ~doc:"Refine MSA after the final Reconciliation step (Default:false)"*)
  +> flag "--mpast"           (optional float)  ~doc:"FLOAT Minimal percentage of alignment of an Amalgam sequences on its (non amalgam) closest sequence to be kept in the final output"
  +> flag "--debug"           no_arg            ~doc:"Get intermediary files (Default:false)"
  +> flag "--just-parse-input"no_arg            ~doc:"Parse input and exit. Recommended to check all input files. (Default:false)"
  +> flag "--html-report"    (optional string)  ~doc:"PATH Logs build events in an HTML report"
  +> flag "--dag-graph"      (optional string)  ~doc:"PATH Write dag graph in an dot file (Can take a lot of time)"
  +> flag "--quiet"           no_arg            ~doc:"Do not report progress.  Default: off"

let command =
  Command.basic
    ~summary:"Amalgam"
    spec
    main

let () = Command.run ~version:"0.1" command

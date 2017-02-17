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

let main sample_sheet outdir species_tree_file alignments_dir seq2sp_dir np memory reconcile refinetree refineali debug () =
  let logger =
    Bistro_logger.tee
      (Bistro_console_logger.create ())
      (Bistro_html_logger.create "report.html")
  in
  let run_reconciliation = Option.value ~default:true reconcile in 
  let refinetree = Option.value ~default:false refinetree in 
  let refineali = Option.value ~default:false refineali in 
  let debug = Option.value ~default:false debug in 
  let np = Option.value ~default:2 np in
  let memory = Option.value ~default:1 memory in
  let configuration = Configuration.load ~sample_sheet ~species_tree_file ~alignments_dir ~seq2sp_dir ~np ~memory ~run_reconciliation ~debug ~refinetree ~refineali ~outdir in
  let amalgam_app = Amalgam.build_app configuration in
  Bistro_app.(
    run ~logger ~np:configuration.Configuration.threads ~mem:(1024 * configuration.Configuration.memory) amalgam_app
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
  +> flag "--reconcile"       (optional bool)   ~doc:"(BOOL if false: Not run final Reconciliation step (Default:true)"
  +> flag "--refinetree"      (optional bool)   ~doc:"(BOOL if true: Refine topology during final Reconciliation step (Default:false)"
  +> flag "--refineali"       (optional bool)   ~doc:"(BOOL if true: Refine MSA after the final Reconciliation step (Default:false)"
  +> flag "--debug"           (optional bool)   ~doc:"(BOOL if true: Get intermediary files (Default:false)"

let command =
  Command.basic
    ~summary:"Amalgam"
    spec
    main

let () = Command.run ~version:"0.1" command

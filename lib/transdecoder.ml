(*
# File: transdecoder.ml
# Created by: Carine Rey
# Created on: July 2016
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

open Core_kernel
open Bistro
open Bistro.Shell_dsl
open Commons



let fasta_template ~fasta ~tmp_fasta =
  let vars = [ "FASTA", dep fasta;
               "FASTA_OK", ident tmp_fasta;
               "NOTEMPTY", string {|`grep -c ">" -m 1 $FASTA`|};
             ]
  in
  let code = {|if [[ "$NOTEMPTY" == "1" ]]; then ln -s $FASTA $FASTA_OK; else echo """>empty
AT
""" > $FASTA_OK; fi|}
  in
  bash_script vars code

let transdecoder
  ?(descr = "")
  ?only_best_orf
  ?only_top_strand
  ?pep_min_length
  ?retain_long_orfs
  ~threads
  ?(memory = 1)
  (fasta:fasta file) : fasta file =
  let descr = if String.is_empty descr then
                  descr
                else
                  ":" ^ descr
  in
  let tmp_fasta = dest // "tmp.fa" in
  Workflow.shell ~descr:("Transdecoder" ^ descr ) ~np:threads ~mem:(Workflow.int (1024 * memory)) [
    mkdir_p dest;
    within_container img (
      and_list [
        mkdir_p tmp;      
        cd dest;
        cmd "bash" [ file_dump (fasta_template ~fasta ~tmp_fasta) ];
        cmd "TransDecoder.LongOrfs" ~img [
          opt "-t" ident tmp_fasta ;
          option (opt "-m" int ) pep_min_length ;
          option (flag string "-S") only_top_strand ;
        ] ;
        cmd "TransDecoder.Predict" ~img [
          opt "-t" ident tmp_fasta ;
          opt "--cpu" ident np ;
          option (flag string "--single_best_orf") only_best_orf ;
          option (opt "--retain_long_orfs_length" int ) retain_long_orfs ;
          string "||";
          string "TransDecoder.Predict";
          opt "-t" ident tmp_fasta ;
          opt "--cpu" ident np ;
          option (flag string "--single_best_orf") only_best_orf ;
          option (opt "--retain_long_orfs_length" int ) retain_long_orfs ;
          string "--no_refine_starts "
        ] ;
        mv (string "*.cds") (string "orfs.cds") ;
      ] ;
    )
  ]
  |> Fn.flip Workflow.select [ "orfs.cds" ]

(*
# File: Generax.ml
# Created by: Carine Rey
# Created on: March 2020
#
#
# Copyright 2020 Carine Rey
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

open Core
open Bistro
open Bistro.Shell_dsl
open Wutils

type generax_configuration = [`generax_configuration] directory

type phylotree


let generax_config ~family ~ali ~link ~tree2 =
  seq ~sep:"\n" [
    seq ~sep:" " [ string "[FAMILIES]"];
    seq ~sep:" " [ string "-"; string family];
    seq ~sep:"=" [ string "starting_gene_tree"; ident tree2];
    seq ~sep:"=" [ string "alignment"; dep ali];
    seq ~sep:"=" [ string "mapping"; dep link];
    seq ~sep:"=" [ string "subst_model"; string "GTR+G"]
  ]

let generax
    ?(descr="")
    ?(memory = 1)
    ~sptreefile
    ~family
    ~threads
    ~link
    ~tree
    (ali :fasta file)
    : phylotree directory =
    let generax_results_prefix = tmp // "generax_output" in
    let output_reconciledtree = dest // (family ^ "_ReconciledTree.nw") in
    let config = tmp // "config.txt" in
    let tree2 = tmp // "tree_without_empty_line.nw" in
    Workflow.shell ~descr:("generax_by_fam" ^ descr) ~version:4 ~np:threads ~mem:(Workflow.int (1024 * memory)) [
    mkdir_p tmp;
    mkdir_p dest;
    within_container caars_img (
      and_list [
        cmd "awk" ~stdout:tree2 [
        string "NF";
        dep tree ]; (*remove empty lines TODO: to be fix earlier*)
        cmd "cat" ~stdout:config [ file_dump (generax_config ~family ~ali ~link ~tree2 ) ];
        cmd "generax" [
          opt "--species-tree" dep sptreefile ;
          opt "--strategy" string "SPR" ; (*Search mode: EVAL does not optimize the tree topology, and just evaluates the likelihood, the DTL rates and the reconciliation of the starting gene trees. SPR performs a tree search (with SPR moves). Default is SPR.*)
          opt "--rec-model" string "UndatedDL" ; (*{UndatedDL, UndatedDTL}; The probabilistic model used to compute the reconciliation likelihood. Default is UndatedDL *)
          opt "--prefix" ident generax_results_prefix;
          opt "--families" ident config;
        ];
        cmd "mv" [
          ident generax_results_prefix // ("results/" ^ family ^ "/geneTree.newick");
          ident output_reconciledtree];
        cmd "python" [
        file_dump (string Scripts.annote_reconciledtree);
        opt "--sp_tree" dep sptreefile;
        opt "--tree" ident output_reconciledtree;
        opt "-sp2seq" dep link;
        opt "--output_prefix" ident (dest // family) ;
        ]
      ]
    )
    ]

(*
# File: profileNJ.ml
# Created by: Carine Rey
# Created on: January 2017
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

open Core
open Bistro.Std
open Bistro.EDSL
open Bistro_bioinfo.Std
open Commons

type phyldog_configuration = [`phyldog_configuration] directory

type phylotree

let script_pre ~tmp_smap ~link =
  let args = [
    "TMP_SMAP", tmp_smap ;
    "LINK", dep link ;
  ]
  in
  bash_script args {|
    cut -f 1 -d ":" $LINK > sp.txt
    cut -f 2 -d ":" $LINK > seq.txt
    paste seq.txt sp.txt > $TMP_SMAP
|}

let script_post ~tree ~tmp_treeout_prefix ~tmp_treeout_prefix2=
  let args = [
    "TREE", dep tree ;
    "TMP_TREEOUT_PREFIX", tmp_treeout_prefix ;
    "TMP_TREEOUT_PREFIX2", tmp_treeout_prefix2 ;
    "DEST", dest ;
  ]
  in
  bash_script args {|
    name=`basename $TREE`
    #get only tree
    cat $TMP_TREEOUT_PREFIX $TMP_TREEOUT_PREFIX2 > trees.txt
    grep -v -h -e ">" trees.txt > only_trees.txt
    echo >> $TREE
    echo >> only_trees.txt
    cat only_trees.txt  $TREE | awk NF > all.tree
    sameroot.py  all.tree all_same_root.tree
    cat  all_same_root.tree | sort | uniq  > $DEST/$name
    |}

let profileNJ
    ?(descr="")
    ~sptreefile
    ~link
    ~tree
    : phylotree directory workflow =

    let threshold_l = [1.0; 0.99; 0.95; 0.9; 0.85] in 
    let tmp_smap = tmp // "smap.txt" in
    let tmp_treeout_prefix = tmp // ("tree_out_pNJ.1.0.tree") in
    let tmp_treeout_prefix2 = tmp // ("tree_out_pNJ.0*.tree") in
    
    let cmd_profileNJ = List.map threshold_l ~f:(fun t ->
    let str_number = Printf.sprintf "%.2f" t in
    let tmp_treeout = tmp // ("tree_out_pNJ." ^ str_number ^ ".tree") in
    cmd "profileNJ" [
      opt "-s" dep sptreefile ;
      opt "-S" ident tmp_smap ;
      opt "-g" dep tree;
      opt "-o" ident tmp_treeout;
      opt "--seuil" float t;
    ]
    ) in

    workflow ~descr:("profileNJ" ^ descr) ~version:4 ~np:1 ~mem:(1024) (List.concat [
    [mkdir_p tmp;
    mkdir_p dest;
    cd tmp;
    (* Preparing profileNJ configuration files*)
    cmd "sh" [ file_dump (script_pre ~tmp_smap ~link) ]];
    cmd_profileNJ;
    [cmd "sh" [ file_dump (script_post ~tree ~tmp_treeout_prefix ~tmp_treeout_prefix2) ]];
    ])

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

let script_pre ~tree ~tmp_treein ~link =
  let args = [
    "TREE", dep tree ;
    "TMP_TREEIN", tmp_treein ;
    "LINK", dep link ;
  ]
  in
  bash_script args {|
    cp $TREE $TMP_TREEIN
    for l in $(cat $LINK); do sp=`echo $l| cut -f 1 -d ":"`; seq=`echo $l| cut -f 2 -d ":"`; sed -i "s/$seq/$seq@$sp/" $TMP_TREEIN; done
|}

let script_post ~tree ~link ~tmp_treeout =
  let args = [
    "TREE", dep tree ;
    "TMP_TREEOUT", tmp_treeout ;
    "LINK", dep link ;
    "DEST", dest ;
  ]
  in
  bash_script args {|
    name=`basename $TREE`
    for l in $(cat $LINK); do sp=`echo $l| cut -f 1 -d ":"`; seq=`echo $l| cut -f 2 -d ":"`; sed -i "s/$seq@$sp/$seq/" $TMP_TREEOUT; done
    grep -v -e ">" $TMP_TREEOUT > $TMP_TREEOUT"_only_tree"
    cat   $TMP_TREEOUT"_only_tree"  $TREE  > $DEST/$name
    |}

let profileNJ
    ?(descr="")
    ~sptreefile
    ~link
    ~tree
    : phylotree directory workflow =

    let threshold_l = [1.0; 0.99; 0.95; 0.9; 0.85] in 
    let tmp_treein = tmp // "tree_in_pNJ.tree" in
    let tmp_treeout = tmp // ("tree_out_pNJ.*.tree") in
    
    let cmd_profileNJ = List.map threshold_l ~f:(fun t ->
    let str_number = Printf.sprintf "%.2f" t in
    let tmp_treeout = tmp // ("tree_out_pNJ." ^ str_number ^ ".tree") in
    cmd "profileNJ" ~env [
      opt "-s" dep sptreefile ;
      opt "-g" ident tmp_treein;
      opt "-o" ident tmp_treeout;
      opt "--seuil" float t;
      opt "--spos" string "postfix";
      opt "--sep" string "@";
    ]
    ) in

    workflow ~descr:("profileNJ" ^ descr) ~version:3 ~np:1 ~mem:(1024) [
    mkdir_p tmp;
    mkdir_p dest;
    cd tmp;
    cmd "sh" [ file_dump (script_pre ~tree ~tmp_treein ~link) ];
    (* Preparing profileNJ configuration files*)
    docker env (
      and_list cmd_profileNJ) ;
    cmd "sh" [ file_dump (script_post ~tree ~link ~tmp_treeout) ];
  ]

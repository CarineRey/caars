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

open Core.Std
open Bistro.Std
open Bistro.EDSL
open Bistro_bioinfo.Std

type phyldog_configuration = [`phyldog_configuration] directory

type phylotree

let profileNJ
    ~sptreefile
    ~threshold
    ~link
    ~tree
    : phylotree directory workflow =

    let tmp_treein = tmp // "tree_in_pNJ.tree" in
    let tmp_treeout = tmp // "tree_in_pNJ.tree" in
    let script_pre = [%bistro {|
    cp {{ dep tree }} {{ ident tmp_treein }}
    for l in $(cat {{ dep link }}); do sp=`echo $l| cut -f 1 -d ":"`; seq=`echo $l| cut -f 2 -d ":"`; sed -i "s/$seq/$seq@$sp/" {{ ident tmp_treein }}; done
    |} ] in
    let script_post = [%bistro {|
    name=`basename {{ dep tree }} `
    for l in $(cat {{ dep link }}); do sp=`echo $l| cut -f 1 -d ":"`; seq=`echo $l| cut -f 2 -d ":"`; sed -i "s/$seq@$sp/$seq/" {{ ident tmp_treeout }} ; done
    grep -v -e ">" {{ ident tmp_treeout }} > {{ ident dest }}/$name
    |} ]
    in

    workflow ~descr:"profileNJ" ~version:3 ~np:1 ~mem:(1024) [
    mkdir_p tmp;
    mkdir_p dest;
    cd tmp;
    (* Preparing profileNJ configuration files*)
    cmd "sh" [ file_dump script_pre ];
    cmd "profileNJ" [
              opt "-s" string sptreefile ;
              opt "-g" ident tmp_treein;
              opt "-o" ident tmp_treeout;
              opt "--seuil" float threshold;
              opt "--spos" string "postfix";
              opt "--sep" string "@";
              ];
    cmd "sh" [ file_dump script_post ];
    ]
    
    

(*
# File: Apytram.ml
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
open Bistro.EDSL_sh
open Bistro_bioinfo.Std
open Commons

let string_of_db_type = function
  | Left F -> "F"
  | Left R -> "R"
  | Left US -> "single"
  | Right RF -> "RF"
  | Right FR -> "FR"
  | Right UP -> "paired"


let apytram 
    ?fastq (* prendre en compte des listes *)
    ?fasta
    ?query
    ?i
    ?no_best_file
    ?only_best_file
    ?evalue
    ?id
    ?fid
    ?mal
    ?fmal
    ?len
    ?flen
    ?required_coverage
    ?stats
    ?plot
    ?plot_ali
    ?(threads = 1)
    ?(memory = 1) 
    ?time_max
    db_type 
    db_blast  : fasta workflow =
    
    let memory = match memory with
      | 0 -> 1
      | _ -> memory
      in 
    workflow ~descr:"apytram.py" ~np:threads ~mem:(memory * 1024) [ (* add memory*)
    cmd "apytram.py" [
        option (opt "-fq" string ) fastq ; (* prendre en compte des listes *)
        option (opt "-fa" string ) fasta ;
        option (opt "-q" dep ) query ;
        option (opt "-i" int ) i ;
        option (opt "-id" float ) id ;
        option (opt "-fid" float ) fid ;
        option (opt "-mal" float ) mal ;
        option (opt "-fmal" float ) fmal ;
        option (opt "-len" float ) len ;
        option (opt "-flen" float ) flen ;
        option (opt "-required_coverage" float ) required_coverage ;
        option (flag string "--stats") stats ;
        option (flag string "--plot") plot ;
        option (flag string "--plot_ali") plot_ali ;
        (*option (opt "-memory" ident) memory ;*)
        (*opt "-memory" int 2 ;*)
        opt "-memory" ident (seq [ string "$((" ; mem ; string " / 1024))" ]) ;
        opt "-threads" ident np ;
        opt "-d" (fun blast_db -> seq [dep db_blast ; string "/db"]) db_blast;
        opt "-dt" string (string_of_db_type db_type);
        opt "-out" seq [ident dest ; string "/apytram"];
        opt "-log" seq [ident dest ; string "/apytram.log"];
        opt "-tmp" ident  tmp ;
        ]
    ]

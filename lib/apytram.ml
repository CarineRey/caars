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
open Bistro.EDSL
open Bistro_bioinfo.Std
open Commons
open Configuration

let string_of_db_type = function
  | Left F -> "F"
  | Left R -> "R"
  | Left US -> "single"
  | Right RF -> "RF"
  | Right FR -> "FR"
  | Right UP -> "paired"

type output

let apytram
    ?fastq
    ?fasta
    ?query
    ?i
    ?no_best_file
    ?only_best_file
    ?write_even_empty
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
    db_blast : output directory workflow =

    let memory = match memory with
      | 0 -> 1
      | _ -> memory
      in

    workflow  ~version:2 ~descr:"apytram.py" ~np:threads ~mem:(memory * 1024) [
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
        option (flag string "--no_best_file") no_best_file ;
        option (flag string "--write_even_empty") write_even_empty ;
        option (flag string "--only_best_file") only_best_file ;
        (*option (opt "-memory" ident) memory ;*)
        (*opt "-memory" int 2 ;*)
        opt "-memory" ident (seq [ string "$((" ; mem ; string " / 1024))" ]) ;
        opt "-threads" ident np ;
        opt "-d" (fun blast_db -> seq [dep blast_db ; string "/db"]) db_blast;
        opt "-dt" string (string_of_db_type db_type);
        opt "-out" seq [ident dest ; string "/apytram"];
        opt "-log" seq [ident dest ; string "/apytram.log"];
        opt "-tmp" ident  ( tmp // "apytram_tmp" ) ;
        (*flag string "--keep_tmp" true;
        opt "-tmp" ident  ( dest // "apytram_tmp" ) ;*)
        ]
    ]

let apytram_multi_species
    ?i
    ?no_best_file
    ?only_best_file
    ?out_by_species
    ?write_even_empty
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
    ~query
    ~fam
    db_blasts : output directory workflow =

    let memory = match memory with
      | 0 -> 1
      | _ -> memory
      in

    let formated_db_blasts =
      List.map db_blasts ~f:(fun (s, w) ->
    seq [dep w ; string "/db:"; string s.id]
      )
    in
    let db_types =
      List.map db_blasts ~f:(fun (s, _) ->
    seq ~sep:":" [string (string_of_db_type (sample_fastq_orientation s.sample_fastq)); string s.id]
      )
    in


    workflow  ~version:4 ~descr:("apytram.py " ^ fam ^ " ")~np:threads ~mem:(memory * 1024) [
    cmd "apytram.py" [
        opt "-q" seq [dep query ; string ":"; string fam] ;
        option (opt "-i" int ) i ;
        option (opt "-id" float ) id ;
        option (opt "-fid" float ) fid ;
        option (opt "-mal" float ) mal ;
        option (opt "-fmal" float ) fmal ;
        option (opt "-len" float ) len ;
        option (opt "-flen" float ) flen ;
        option (opt "-required_coverage" float ) required_coverage ;
        option (flag string "--stats") stats ;
        option (flag string "--no_best_file") no_best_file ;
        option (flag string "--write_even_empty") write_even_empty ;
        option (flag string "--only_best_file") only_best_file ;
        option (flag string "--out_by_species") out_by_species ;
        opt "-memory" ident (seq [ string "$((" ; mem ; string " / 1024))" ]) ;
        opt "-threads" ident np ;
        opt "-d" ident (seq ~sep:"," formated_db_blasts) ;
        opt "-dt" ident (seq ~sep:"," db_types) ;
        opt "-out" seq [ident dest ; string "/apytram"] ;
        opt "-log" seq [ident dest ; string "/apytram.log"] ;
        opt "-tmp" ident  ( tmp // "apytram_tmp" ) ;
        (*flag string "--keep_tmp" true;
        opt "-tmp" ident  ( dest // "apytram_tmp" ) ;*)
        ]
    ]


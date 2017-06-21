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

open Core
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

type apytram_output

let apytram_multi_species
    ?(descr="")
    ?i
    ?evalue
    ?no_best_file
    ?only_best_file
    ?out_by_species
    ?write_even_empty
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
    (compressed_reads_dbs : compressed_read_db list) : apytram_output directory workflow =

    let memory = match memory with
      | 0 -> 1
      | _ -> memory
      in

    let formated_db_blasts =
      List.map compressed_reads_dbs ~f:(fun db ->
    seq [dep db.cluster_rep_blast_db ; string "/db:"; string db.s.id]
      )
    in
    let db_types =
      List.map compressed_reads_dbs ~f:(fun {s}->
    seq ~sep:":" [string (string_of_db_type (sample_file_orientation s.sample_file)); string s.id]
      )
    in
    let formated_fasta =
      List.map compressed_reads_dbs ~f:(fun db ->
    seq [dep db.concat_fasta ; string ":"; string db.s.id]
      )
    in
    let formated_fastaidx =
      List.map compressed_reads_dbs ~f:(fun db ->
    seq [dep db.index_concat_fasta ; string ":"; string db.s.id]
      )
    in
    let formated_cluster =
      List.map compressed_reads_dbs ~f:(fun db ->
    seq [dep db.reformated_cluster ; string ":"; string db.s.id]
      )
    in
    let formated_clusteridx =
      List.map compressed_reads_dbs ~f:(fun db ->
    seq [dep db.index_cluster ; string ":"; string db.s.id]
      )
    in


    workflow  ~version:5 ~descr:("apytram.py" ^ descr) ~np:threads ~mem:(memory * 1024) [
    cmd "apytram.py" [
        opt "-q" seq [dep query ; string ":"; string fam] ;
        option (opt "-i" int ) i ;
        option (opt "-e" float ) evalue;
        option (opt "-id" float ) id ;
        option (opt "-fid" float ) fid ;
        option (opt "-mal" float ) mal ;
        option (opt "-fmal" float ) fmal ;
        option (opt "-len" float ) len ;
        option (opt "-flen" float ) flen ;
        option (opt "-required_coverage" float ) required_coverage ;
        option (opt "-time_max" int ) time_max ;
        option (flag string "--stats") stats ;
        option (flag string "--no_best_file") no_best_file ;
        option (flag string "--write_even_empty") write_even_empty ;
        option (flag string "--only_best_file") only_best_file ;
        option (flag string "--out_by_species") out_by_species ;
        opt "-memory" ident (seq [ string "$((" ; mem ; string " / 1024))" ]) ;
        opt "-threads" ident np ;
        opt "-d" ident (seq ~sep:"," formated_db_blasts) ;
        opt "-dt" ident (seq ~sep:"," db_types) ;
        opt "-fa" ident (seq ~sep:"," formated_fasta) ;
        opt "-idx" ident (seq ~sep:"," formated_fastaidx) ;
        flag string "--UseIndex" true;
        opt "-clstr" ident (seq ~sep:"," formated_cluster) ;
        opt "-clstridx" ident (seq ~sep:"," formated_clusteridx) ;
        opt "-out" seq [ident dest ; string "/apytram"] ;
        opt "-log" seq [ident dest ; string "/apytram.log"] ;
        opt "-tmp" ident  ( tmp // "apytram_tmp" ) ;
        (*flag string "--keep_tmp" true;
        opt "-tmp" ident  ( dest // "apytram_tmp" ) ;*)
        ]
    ]


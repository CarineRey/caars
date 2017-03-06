(*
# File: trinity.ml
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


type assembly_stats

let ( % ) f g = fun x -> g (f x)



let single_stranded_or_unstranded = function
  | F -> string "--SS_lib_type F"
  | R -> string "--SS_lib_type R"
  | US -> string ""

let paired_stranded_or_unstranded = function
  | RF -> string "--SS_lib_type RF"
  | FR -> string "--SS_lib_type RR"
  | UP -> string ""

let config_trinity_fasta_paired_or_single = function
  | Fasta_Single_end (w, o ) ->
        seq ~sep: " " [ string "--single" ; dep w ; single_stranded_or_unstranded o ]
  | Fasta_Paired_end (lw, rw , o) ->
       seq ~sep: " " [ string "--left" ; dep lw ; string "--right" ; dep rw ; paired_stranded_or_unstranded o]

let trinity_fasta
    ?(descr = "")
    ?full_cleanup
    ?no_normalization
    ~threads
    ?(memory = 1)
    (sample_fasta : fasta workflow sample_fasta)
    : fasta workflow =
    let descr = if descr = "" then
                  descr
                else
                  ":" ^ descr ^ " "
    in
    workflow ~descr:("Trinity" ^ descr) ~np:threads ~mem:(1024 * memory) [
        mkdir_p dest;
        cmd "Trinity" [
            string "--no_version_check";
            opt "--max_memory" ident (seq [ string "$((" ; mem ; string " / 1024))G" ]) ;
            opt "--CPU" ident np ;
            option (flag string "--full_cleanup") full_cleanup ;
            option (flag string "--no_normalize_reads") no_normalization ;
            config_trinity_fasta_paired_or_single sample_fasta;
            string "--seqType fa" ;
            opt "--output" seq [ ident dest ; string "/trinity"] ;
        ];
        cmd "sed" [
                     string "-re";
                     string {|"s/(>[_a-zA-Z0-9]*)( len=[0-9]* path=.*)/\1/"|};
                     string "-i";
                     seq [ident dest; string "/trinity.Trinity.fasta";];
        ];
    ]
    / selector [ "trinity.Trinity.fasta" ]

let config_fastq_paired_or_single = function
  | Fastq_Single_end (w, _ ) ->
        seq ~sep: " " [ string "--single" ; dep w ]
  | Fastq_Paired_end (lw, rw , _) ->
       seq ~sep: " " [ string " --pairs_together" ;  string "--left" ; dep lw; string "--right" ; dep rw; string "--pairs_together --PARALLEL_STATS" ]

let config_output_fastq_paired_or_single = function
  | Fastq_Single_end (w, _ ) ->
       seq ~sep: "\n" [  string "singlelink=`readlink single.norm.fq`";
                         seq ~sep: " " [string "mv $singlelink";  dest // "single.norm.fq" ];
                      ]
  | Fastq_Paired_end (lw, rw , _) ->
       seq ~sep: "\n" [  string "leftlink=`readlink left.norm.fq`";
                         seq ~sep: " " [string "mv $leftlink";  dest // "left.norm.fq" ];
                         string "rightlink=`readlink right.norm.fq`";
                         seq ~sep: " " [string "mv $rightlink";  dest // "right.norm.fq" ];
                      ]

let fastq_read_normalization
    max_cov
    ~threads
    ?(memory = 1)
    (fastq : _ fastq workflow sample_fastq)
    : _ fastq directory workflow =
  let script = [%bistro{|
    TRINITY_PATH=`which Trinity`
    TRINTIY_DIR_PATH=`dirname $TRINITY_PATH`
    READ_NORMALISATION_PATH=$TRINTIY_DIR_PATH/util/insilico_read_normalization.pl
    $READ_NORMALISATION_PATH  {{ config_fastq_paired_or_single fastq }} --seqType "fq" --JM {{ seq [ string "$((" ; mem ; string " / 1024))" ]}}G --max_cov {{ int max_cov }} --CPU {{ ident np }} --output {{ ident tmp }}

    {{ config_output_fastq_paired_or_single fastq }}

    |}]
  in
  workflow ~descr:"fastq_read_normalization" ~version:2 ~np:threads ~mem:(1024 * memory) [
    mkdir_p dest;
    mkdir_p tmp;
    cd tmp;
    cmd "sh" [ file_dump script ]
    ]
let config_fasta_paired_or_single = function
  | Fasta_Single_end (w, _ ) ->
        seq ~sep: " " [ string "--single" ; dep w ]
  | Fasta_Paired_end (lw, rw , _) ->
       seq ~sep: " " [ string " --pairs_together" ;  string "--left" ; dep lw; string "--right" ; dep rw; string "--pairs_together --PARALLEL_STATS" ]

let config_output_fasta_paired_or_single = function
  | Fasta_Single_end (w, _ ) ->
       seq ~sep: "\n" [  string "singlelink=`readlink single.norm.fa`";
                         seq ~sep: " " [string "mv $singlelink";  dest // "single.norm.fa" ];
                      ]
  | Fasta_Paired_end (lw, rw , _) ->
       seq ~sep: "\n" [  string "leftlink=`readlink left.norm.fa`";
                         seq ~sep: " " [string "mv $leftlink";  dest // "left.norm.fa" ];
                         string "rightlink=`readlink right.norm.fa`";
                         seq ~sep: " " [string "mv $rightlink";  dest // "right.norm.fa" ];
                      ]

let fasta_read_normalization
    ?(descr = "")
    max_cov
    ~threads
    ?(memory = 1)
    (fasta : fasta workflow sample_fasta)
    : fasta directory workflow =
  let descr = if descr = "" then
                  descr
                else
                  ":" ^ descr ^ " "
  in
  let script = [%bistro{|
    TRINITY_PATH=`which Trinity`
    TRINTIY_DIR_PATH=`dirname $TRINITY_PATH`
    READ_NORMALISATION_PATH=$TRINTIY_DIR_PATH/util/insilico_read_normalization.pl
    $READ_NORMALISATION_PATH  {{ config_fasta_paired_or_single fasta }} --seqType "fa" --JM {{ seq [ string "$((" ; mem ; string " / 1024))" ]}}G --max_cov {{ int max_cov }} --CPU {{ ident np }} --output {{ ident tmp }}

    {{ config_output_fasta_paired_or_single fasta }}

    |}]
  in
  workflow ~descr:("fasta_read_normalization" ^ descr) ~version:2 ~np:threads ~mem:(1024 * memory) [
    mkdir_p dest;
    mkdir_p tmp ;
    cd tmp;
    cmd "sh" [ file_dump script ];
    ]


let fastool ?(descr="") ~dep_input (fastq : _ fastq workflow) :  fasta workflow =
  let descr = if descr = "" then
                  descr
                else
                  ":" ^ descr ^ " "
  in
  let script = [%bistro {|
    TRINITY_PATH=`which Trinity`
    TRINTIY_DIR_PATH=`dirname $TRINITY_PATH`
    FASTOOL_PATH=$TRINTIY_DIR_PATH/trinity-plugins/fastool/fastool
    $FASTOOL_PATH --illumina-trinity --to-fasta  {{ dep fastq }} > {{ ident dest }}
    |} ]
  in
  workflow ~descr:("fastq2fasta" ^ descr) ~np:1 [
    cmd "ls" [ dep dep_input ];
    cmd "sh" [ file_dump script ];
  ]


let assembly_stats ?(descr="") (fasta:fasta workflow) : assembly_stats workflow =
   let descr = if descr = "" then
                  descr
                else
                  ":" ^ descr ^ " "
   in
   let script = [%bistro {|
    TRINITY_PATH=`which Trinity`
    TRINTIY_DIR_PATH=`dirname $TRINITY_PATH`
    TRINITYSTATS_PATH=$TRINTIY_DIR_PATH/util/TrinityStats.pl
    FASTA={{ dep fasta }}
    if [ -s $FASTA ]
    then
    $TRINITYSTATS_PATH {{ dep fasta }} > {{ ident dest }}
    else
    echo "Empty file" > {{ ident dest }}
    fi
    |} ]
  in
  workflow ~descr:("assembly_stats_trinity" ^ descr) ~np:1 [
    cmd "sh" [ file_dump script ];
  ]

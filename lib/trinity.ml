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


let ( % ) f g = fun x -> g (f x)


let trinity_fasta
    ?full_cleanup
    ~threads
    ?(memory = 1)
    ~is_paired
    (fasta: fasta workflow) : fasta workflow =
    workflow  ~np:threads ~mem:(1024 * memory) [
        mkdir_p dest;
        cmd "Trinity" [
            opt "--max_memory" ident (seq [ string "$((" ; mem ; string " / 1024))G" ]) ;
                    opt "--CPU" ident np ;
                    option (flag string "--full_cleanup") full_cleanup ;
                    opt "-single" dep fasta;
                    flag string "--run_as_paired" is_paired;
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
    fastq : _ fastq directory workflow =
    workflow ~version:2 ~np:threads ~mem:(1024 * memory) [
    mkdir_p dest;
    mkdir_p tmp;
    cd tmp;
    script "sh" [%bistro{|
    TRINITY_PATH=`which Trinity`
    TRINTIY_DIR_PATH=`dirname $TRINITY_PATH`
    READ_NORMALISATION_PATH=$TRINTIY_DIR_PATH/util/insilico_read_normalization.pl
    $READ_NORMALISATION_PATH  {{ config_fastq_paired_or_single fastq }} --seqType "fq" --JM {{ seq [ string "$((" ; mem ; string " / 1024))" ]}}G --max_cov {{ int max_cov }} --CPU {{ ident np }} --output {{ ident tmp }}

    {{ config_output_fastq_paired_or_single fastq }}

    |}]
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
    max_cov
    ~threads
    ?(memory = 1)
    (fasta : fasta workflow sample_fasta)
    : fasta directory workflow =
    workflow ~version:2 ~np:threads ~mem:(1024 * memory) [
    mkdir_p dest;
    mkdir_p tmp ;
    cd tmp;
    script "sh" [%bistro{|
    TRINITY_PATH=`which Trinity`
    TRINTIY_DIR_PATH=`dirname $TRINITY_PATH`
    READ_NORMALISATION_PATH=$TRINTIY_DIR_PATH/util/insilico_read_normalization.pl
    $READ_NORMALISATION_PATH  {{ config_fasta_paired_or_single fasta }} --seqType "fa" --JM {{ seq [ string "$((" ; mem ; string " / 1024))" ]}}G --max_cov {{ int max_cov }} --CPU {{ ident np }} --output {{ ident tmp }}

    {{ config_output_fasta_paired_or_single fasta }}

    |}]
    ]


let fastool (fastq : _ fastq workflow) :  fasta workflow =
    workflow ~np:1 [script "sh" [%bistro {|
    TRINITY_PATH=`which Trinity`
    TRINTIY_DIR_PATH=`dirname $TRINITY_PATH`
    FASTOOL_PATH=$TRINTIY_DIR_PATH/trinity-plugins/fastool/fastool
    $FASTOOL_PATH --illumina-trinity --to-fasta  {{ dep fastq }} > {{ ident dest }}
    |} ]
    ]

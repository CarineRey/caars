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

open Core
open Bistro
open Bistro.Shell_dsl
open Commons


type assembly_stats

let ( % ) f g = fun x -> g (f x)



let single_stranded_or_unstranded = function
  | F -> string "--SS_lib_type F"
  | R -> string "--SS_lib_type R"
  | US -> string ""

let paired_stranded_or_unstranded = function
  | RF -> string "--SS_lib_type RF"
  | FR -> string "--SS_lib_type FR"
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
    (sample_fasta : fasta pworkflow sample_fasta)
    : fasta pworkflow =
    let descr = if String.is_empty descr then
                  descr
                else
                  ":" ^ descr
    in
    Workflow.shell ~descr:("Trinity" ^ descr) ~np:threads ~mem:(Workflow.int (1024 * memory)) [
        mkdir_p dest;
        cmd "Trinity" ~img [
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
    |> Fn.flip Workflow.select [ "trinity.Trinity.fasta" ]

let config_fastq_paired_or_single = function
  | Fastq_Single_end (w, _ ) ->
        seq ~sep: " " [ string "--single" ; dep w ]
  | Fastq_Paired_end (lw, rw , _) ->
       seq ~sep: " " [ string " --pairs_together" ;  string "--left" ; dep lw; string "--right" ; dep rw; string "--pairs_together --PARALLEL_STATS" ]

let config_output_fastq_paired_or_single = function
  | Fastq_Single_end _ ->
       seq ~sep: "\n" [  string "singlelink=`readlink single.norm.fq`";
                         seq ~sep: " " [string "mv $singlelink";  dest // "single.norm.fq" ];
                      ]
  | Fastq_Paired_end _ ->
       seq ~sep: "\n" [  string "leftlink=`readlink left.norm.fq`";
                         seq ~sep: " " [string "mv $leftlink";  dest // "left.norm.fq" ];
                         string "rightlink=`readlink right.norm.fq`";
                         seq ~sep: " " [string "mv $rightlink";  dest // "right.norm.fq" ];
                      ]

(* let fastq_read_normalization *)
(*     max_cov *)
(*     ~threads *)
(*     ?(memory = 1) *)
(*     (fastq : _ fastq workflow sample_fastq) *)
(*     : _ fastq directory workflow = *)
(*   let script = [%bistro{| *)
(*     TRINITY_PATH=`which Trinity` *)
(*     TRINTIY_DIR_PATH=`dirname $TRINITY_PATH` *)
(*     READ_NORMALISATION_PATH=$TRINTIY_DIR_PATH/util/insilico_read_normalization.pl *)
(*     $READ_NORMALISATION_PATH  {{ config_fastq_paired_or_single fastq }} --seqType "fq" --JM {{ seq [ string "$((" ; mem ; string " / 1024))" ]}}G --max_cov {{ int max_cov }} --CPU {{ ident np }} --output {{ ident tmp }} *)

(*     {{ config_output_fastq_paired_or_single fastq }} *)

(*     |}] *)
(*   in *)
(*   workflow ~descr:"fastq_read_normalization" ~version:2 ~np:threads ~mem:(1024 * memory) [ *)
(*     mkdir_p dest; *)
(*     mkdir_p tmp; *)
(*     cd tmp; *)
(*     cmd "sh" [ file_dump script ] *)
(*     ] *)

let config_fasta_paired_or_single = function
  | Fasta_Single_end (w, _ ) ->
        seq ~sep: " " [ string "--single" ; dep w ]
  | Fasta_Paired_end (lw, rw , _) ->
       seq ~sep: " " [ string " --pairs_together" ;  string "--left" ; dep lw; string "--right" ; dep rw; string "--pairs_together --PARALLEL_STATS" ]


let config_output_fasta_paired_or_single = function
  | Fasta_Single_end _ ->
       seq ~sep: "\n" [  string "singlelink=`readlink single.norm.fa`";
                         seq ~sep: " " [string "mv $singlelink";  dest // "single.norm.fa" ];
                      ]
  | Fasta_Paired_end _ ->
       seq ~sep: "\n" [  string "leftlink=`readlink left.norm.fa`";
                         seq ~sep: " " [string "mv $leftlink";  dest // "left.norm.fa" ];
                         string "rightlink=`readlink right.norm.fa`";
                         seq ~sep: " " [string "mv $rightlink";  dest // "right.norm.fa" ];
                      ]

(* let fasta_read_normalization *)
(*     ?(descr = "") *)
(*     max_cov *)
(*     ~threads *)
(*     ?(memory = 1) *)
(*     (fasta : fasta workflow sample_fasta) *)
(*     : fasta directory workflow = *)
(*   let descr = if descr = "" then *)
(*                   descr *)
(*                 else *)
(*                   ":" ^ descr ^ " " *)
(*   in *)
(*   let memory = int_of_float( float_of_int memory *. 0.75) in *)
(*   let script = [%bistro{| *)
(*     TRINITY_PATH=`which Trinity` *)
(*     TRINTIY_DIR_PATH=`dirname $TRINITY_PATH` *)
(*     READ_NORMALISATION_PATH=$TRINTIY_DIR_PATH/util/insilico_read_normalization.pl *)
(*     $READ_NORMALISATION_PATH  {{ config_fasta_paired_or_single fasta }} --seqType "fa" --JM {{ seq [ string "$((" ; mem ; string " / 1024))" ]}}G --max_cov {{ int max_cov }} --CPU {{ ident np }} --output {{ ident tmp }} *)

(*     {{ config_output_fasta_paired_or_single fasta }} *)

(*     |}] *)
(*   in *)
(*   workflow ~descr:("fasta_read_normalization_(custom)" ^ descr) ~version:2 ~np:threads ~mem:(1024 * memory) [ *)
(*     mkdir_p dest; *)
(*     mkdir_p tmp ; *)
(*     cd tmp; *)
(*     cmd "sh" [ file_dump script ]; *)
(*     ] *)

let fasta_read_normalization_script ~fasta ~max_cov ~given_mem=
  let vars = [
    "TRINITY_PATH", string "`which Trinity`" ;
    "TRINTIY_DIR_PATH", string "`dirname $TRINITY_PATH`" ;
    "PERL5LIB", string "$TRINTIY_DIR_PATH/PerlLib:$PERL5LIB" ;
    "FA", config_fasta_paired_or_single fasta ;
    "MEM", seq [ string "$((" ; int given_mem ; string " / 1024))" ] ;
    "MAX_COV", int max_cov ;
    "NP", np ;
    "TMP", tmp ;
  ]
  in
  bash_script vars {|
    trinity.insilico_read_normalization.pl $FA --seqType "fa" --JM ${MEM}G --max_cov $MAX_COV --CPU $NP --output $TMP --trinity_dir $TRINTIY_DIR_PATH
    |}



let fasta_read_normalization_get_output ~fasta ~dest=
  let (vars, code) = match fasta with
    | Fasta_Single_end _ -> (["DEST", dest;
                              "SINGLELINK", string "`readlink single.norm.fa`"],
                             {| mv $SINGLELINK $DEST/"single.norm.fa"|})
    | Fasta_Paired_end _ -> (["DEST", dest;
                              "LEFTLINK", string "`readlink left.norm.fa`";
                              "RIGHTLINK", string "`readlink right.norm.fa`"],
                             {|echo $LEFTLINK ; mv $LEFTLINK $DEST/"left.norm.fa"; mv $RIGHTLINK $DEST/"right.norm.fa"|})
  in
  bash_script vars code

let fasta_read_normalization_2
    ?(descr = "")
    max_cov
    ~threads
    ?(memory = 1)
    ?(max_memory = 1)
    (fasta : fasta pworkflow sample_fasta)
    : fasta dworkflow =
  let descr = if String.is_empty "" then
                  descr
                else
                  ":" ^ descr
  in

  let bistro_memory = if max_memory > 2
                      then
                         Stdlib.(min max_memory (int_of_float( float_of_int memory *. 2.)))
                      else
                         1
                      in
  let given_mem =    if bistro_memory > 2
                      then
                         Stdlib.(int_of_float( float_of_int bistro_memory /. 2.))
                      else
                         1
                      in
  (* reserve more memory by bistro than given to normalization tools*)
  Workflow.shell ~descr:("fasta_read_normalization_(custom)" ^ descr) ~version:2 ~np:threads ~mem:(Workflow.int (1024 * bistro_memory)) [
    mkdir_p dest;
    mkdir_p tmp ;
    within_container img (
      and_list [
        cd tmp ;
        cmd "sh" [ file_dump (fasta_read_normalization_script ~fasta ~max_cov ~given_mem) ];
        cmd "sh" [ file_dump (fasta_read_normalization_get_output ~fasta ~dest) ];
      ]
    )
  ]

let fastq2fasta ?(descr="") ?(dep_input=None) (fastq : #fastq pworkflow) :  fasta pworkflow =
    let (check_input, w_input) = match dep_input with
                        | Some w -> (true, dep w)
                        | None -> (false, ident dest)
                        in
    let descr = if String.is_empty descr then descr else ":" ^ descr in
    let script =
        let vars = [
        "CHECK", flag seq [string "ls "; w_input] check_input  ;
        "FQ", dep fastq ;
        "DEST", dest ;
        ]
        in
        bash_script vars {|
        $CHECK
        seqtk seq -A $FQ > $DEST
        |}
    in
    Workflow.shell ~descr:("fastq2fasta" ^ descr) ~np:1 [
        cmd "sh" ~img [ file_dump script ];
    ]

let assembly_stats ?(descr="") (fasta:fasta pworkflow) : assembly_stats pworkflow =
   let descr = if String.is_empty descr then descr else ":" ^ descr ^ " " in
   let script =
     let vars = [
       "TRINITY_PATH", string "`which Trinity`" ;
       "TRINTIY_DIR_PATH", string "`dirname $TRINITY_PATH`" ;
       "TRINITYSTATS_PATH", string "$TRINTIY_DIR_PATH/util/TrinityStats.pl" ;
       "FASTA" , dep fasta ;
       "DEST", dest ;
     ]
     in
     bash_script vars {|
    if [ -s $FASTA ]
    then
    $TRINITYSTATS_PATH $FASTA > $DEST
    else
    echo "Empty file" > $DEST
    fi
    |}
  in
  Workflow.shell ~descr:("assembly_stats_trinity" ^ descr) ~np:1 [
    cmd "sh" ~img [ file_dump script ];
  ]

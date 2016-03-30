open Core.Std
open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std



(* memory bound correspond to storing a human index in memory, following bowtie manual *)
let apytram 
    ?fastq (* prendre en compte des listes *)
    ?fasta
    ?query
    ?tmp
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
    ?threads
    ?memory
    ?time_max
    db_type 
    db_blast  : fasta workflow =
    
    workflow ~descr:"apytram.py" ~mem:(3 * 1024) ?np:threads [ (* add memory*)
    cmd "apytram.py" [
        option (opt "-fq" string ) fastq ; (* prendre en compte des listes *)
        option (opt "-fa" string ) fasta ;
        option (opt "-q" string ) query ;
        option (opt "-tmp" string ) tmp ;
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
        option (opt "-memory" int) memory ;
        option (opt "-threads" int) threads ;
        opt "-d" (fun blast_db -> seq [string db_blast ; string "/db"]) db_blast;
        opt "-dt" string db_type ; (* Check {single,paired,FR,RF,F,R}*)
        opt "-out" seq [ident dest ; string "/apytram"];
        opt "-log" seq [ident dest ; string "/apytram.log"]
        ]
    ]
    
(**
-d DATABASE -dt {single,paired,FR,RF,F,R}
                  [-fa [FASTA [FASTA ...]]] [-fq [FASTQ [FASTQ ...]]]
                  [-q QUERY] [-pep QUERY_PEP] [-i ITERATION_MAX]
                  [-i_start ITERATION_START] [-out OUTPUT_PREFIX] [-log LOG]
                  [-tmp TMP] [--keep_iterations] [--no_best_file]
                  [--only_best_file] [--stats] [--plot] [--plot_ali]
                  [-e EVALUE] [-id MIN_ID] [-mal MIN_ALI_LEN] [-len MIN_LEN]
                  [-required_coverage REQUIRED_COVERAGE] [--finish_all_iter]
                  [-flen FINAL_MIN_LEN] [-fid FINAL_MIN_ID]
                  [-fmal FINAL_MIN_ALI_LEN] [-threads THREADS]
                  [-memory MEMORY] [-time_max TIME_MAX]




let qual_option (type s) x = match (x : s Fastq.format) with
  | Fastq.Solexa  -> "--solexa-quals"
  | Fastq.Sanger -> "--phred33-quals"
  | Fastq. Phred64 -> "--phred64-quals"

let flag_of_preset mode preset =
  let flag = match preset with
    | `very_fast -> "--very-fast"
    | `fast -> "--fast"
    | `sensitive -> "--sensitive"
    | `very_sensitive -> "--very-sensitive"
  in
  if mode = `local then flag ^ "-local" else flag

let flag_of_mode = function
  | `end_to_end -> "--end-to-end"
  | `local -> "--local"

let flag_of_db_type = function
  | `single -> "single"
  | `paired -> "paired"
  | `FR -> "FR"
  | `RF -> "RF"
  | `F -> "F"
  | `R -> "-R"
  
let bowtie2
    ?skip ?qupto ?trim5 ?trim3 ?preset
    ?_N ?_L ?ignore_quals ?(mode = `end_to_end)
    ?a ?k ?_D ?_R ?minins ?maxins ?orientation
    ?no_mixed ?no_discordant ?dovetail ?no_contain ?no_overlap
    ?threads ?seed
    ?fastq_format index fqs =

  let args = match fqs with
    | `single_end fqs ->
      opt "-U" (list dep ~sep:",") fqs
    | `paired_end (fqs1, fqs2) ->
      seq [
        opt "-1" (list dep ~sep:",") fqs1 ;
        string " " ;
        opt "-2" (list dep ~sep:",") fqs2
      ]
  in
  workflow ~descr:"bowtie2" ~mem:(3 * 1024) ?np:threads ~pkgs:[package] [
    cmd "bowtie2" [
      option (opt "--skip" int) skip ;
      option (opt "--qupto" int) qupto ;
      option (opt "--trim5" int) trim5 ;
      option (opt "--trim3" int) trim3 ;
      option ((flag_of_preset mode) % string) preset ;
      option (opt "-N" int) _N ;
      option (opt "-L" int) _L ;
      option (flag string "--ignore-quals") ignore_quals ;
      (flag_of_mode % string) mode ;
      option (flag string "-a") a ;
      option (opt "-k" int) k ;
      option (opt "-D" int) _D ;
      option (opt "-R" int) _R ;
      option (opt "--minins" int) minins ;
      option (opt "--maxins" int) maxins ;
      option (flag_of_orientation % string) orientation ;
      option (flag string "--no-mixed") no_mixed  ;
      option (flag string "--no-discordant") no_overlap  ;
      option (flag string "--dovetail") dovetail ;
      option (flag string "--no-contain") no_contain ;
      option (flag string "--no-overlap") no_overlap ;
      option (opt "--threads" int) threads ;
      option (opt "--seed" int) seed ;
      option (opt "-q" (qual_option % string)) fastq_format ;
      opt "-x" (fun index -> seq [dep index ; string "/index"]) index ;
      args ;
      opt "-S" ident dest ;
    ]
  ]

**)

open Core.Std
open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std


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
    ?threads
    ?memory
    ?time_max
    db_type 
    db_blast  : fasta workflow =
    
     
    workflow ~descr:"apytram.py"  ~mem:(3 * 1024) ?np:threads [ (* add memory*)
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
        option (opt "-memory" int) memory ;
        option (opt "-threads" int) threads ;
        opt "-d" (fun blast_db -> seq [dep db_blast ; string "/db"]) db_blast;
        opt "-dt" string db_type ; (* Check {single,paired,FR,RF,F,R}*)
        opt "-out" seq [ident dest ; string "/apytram"];
        opt "-log" seq [ident dest ; string "/apytram.log"];
        opt "-tmp" ident  tmp ;
        ]
    ]

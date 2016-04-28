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

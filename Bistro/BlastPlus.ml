open Core_kernel.Std
open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std



type blast_db = [`blast_db] directory
 
let makeblastdb ?parse_seqids ~dbtype  dbname  (fasta : fasta workflow)  : blast_db workflow =
	workflow [
		mkdir_p dest;
		cmd "makeblastdb" [ option (flag string "-parse_seqids") parse_seqids ;
				    opt "-in" dep fasta;
				    opt "-dbtype" string dbtype ;
				    string "-out" ; seq ~sep:"/" [ dest ;string dbname; string "db" ] ] ;
		 ]
   / selector [ dbname ] 
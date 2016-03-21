open Core_kernel.Std
open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std

let trinity (fastq: _ fastq workflow) : fasta workflow =
	workflow [
		mkdir_p dest;
		cmd "Trinity" [ string "-single" ; dep fastq ; string "--seqType fq --max_memory 1G --full_cleanup " ; 
		string "--output" ; seq ~sep:"/" [ dest ; string "trinity"] ] ;
	]
    / selector [ "trinity.Trinity.fasta" ] 



open Core.Std
open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std

let ( % ) f g = fun x -> g (f x)


let run_as_paired x = match x with
  | "single"  -> ""
  | "paired" -> "--run_as_paired"
  | _ -> ""
  
let trinity_fasta
	?threads
	?full_cleanup
	~memory
	~fasta_type
	(fasta: fasta workflow) : fasta workflow =
	workflow ?np:threads [
		mkdir_p dest;
		cmd "Trinity" [
		    opt "--max_memory" (sprintf "%dG" % string) memory ;
            option (opt "--CPU" int) threads ;
            option (flag string "--full_cleanup") full_cleanup ;
            opt "-single" dep fasta;
            string (run_as_paired fasta_type);
			string "--seqType fa" ;
			opt "--output" seq [ ident dest ; string "/trinity"] ;
			]
	]
    / selector [ "trinity.Trinity.fasta" ]


(*let fastq2fasta (fastq : _ fastq workflow ) : fasta workflow = 
       workflow [
            cmd "awk" ~stdout:dest [ string {|"NR%4==1||NR%4==2"|} ; dep fastq ; string {| | tr @ ">" |}]
    ]
*)
    
let read_normalization seq_type memmory max_cov nb_cpu (fastq : _ fastq workflow) : fasta workflow =
	Bistro.Workflow.make [%sh{|
	TRINITY_PATH=`which Trinity`
	TRINTIY_DIR_PATH=`dirname $TRINITY_PATH`
	READ_NORMALISATION_PATH=$TRINTIY_DIR_PATH/util/insilico_read_normalization.pl
	$READ_NORMALISATION_PATH --single {{ dep fastq }} --seqType {{ string seq_type }} --JM {{ int memmory}}G --max_cov {{ int max_cov}} --CPU {{int nb_cpu}} --output {{ ident dest }} --no_cleanup
	|}]
	/ selector [ "tmp_normalized_reads/single.fa"]
	
	
let fastool (fastq : _ fastq workflow) :  fasta workflow =
	Bistro.Workflow.make [%sh {|
	TRINITY_PATH=`which Trinity`
	TRINTIY_DIR_PATH=`dirname $TRINITY_PATH`
	FASTOOL_PATH=$TRINTIY_DIR_PATH/trinity-plugins/fastool/fastool
	$FASTOOL_PATH --illumina-trinity --to-fasta  {{ dep fastq }} > {{ ident dest }}
	|} ]


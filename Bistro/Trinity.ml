open Core.Std
open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std

let trinity (fastq: _ fastq workflow) : fasta workflow =
	workflow [
		mkdir_p dest;
		cmd "Trinity" [ opt "-single" dep fastq;
						string "--seqType fq --max_memory 1G --full_cleanup " ;
						string "--output" ; seq ~sep:"/" [ dest ; string "trinity"] ] ;
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
	$READ_NORMALISATION_PATH --single {{ dep fastq }} --seqType {{ string seq_type }} --JM {{ string memmory}}G --max_cov {{ string max_cov}} --CPU {{string nb_cpu}} --output {{ ident dest }} --no_cleanup
	|}]
	/ selector [ "tmp_normalized_reads/single.fa"]
	
	
let fastool (fastq : _ fastq workflow) :  fasta workflow =
	Bistro.Workflow.make [%sh {|
	TRINITY_PATH=`which Trinity`
	TRINTIY_DIR_PATH=`dirname $TRINITY_PATH`
	FASTOOL_PATH=$TRINTIY_DIR_PATH/trinity-plugins/fastool/fastool
	$FASTOOL_PATH --illumina-trinity --to-fasta  {{ dep fastq }} > {{ ident dest }}
	|} ]


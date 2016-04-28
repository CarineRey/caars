open Core.Std
open Bistro.Std
open Bistro.EDSL_sh
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
		]
	]
    / selector [ "trinity.Trinity.fasta" ]


(*let fastq2fasta (fastq : _ fastq workflow ) : fasta workflow = 
       workflow [
            cmd "awk" ~stdout:dest [ string {|"NR%4==1||NR%4==2"|} ; dep fastq ; string {| | tr @ ">" |}]
    ]
*)

let config_paired_or_single = function
  | Single_end (w, _ ) ->
        seq ~sep: " " [ string "--single" ; dep w ]
  | Paired_end (lw, rw , _) ->
       seq ~sep: " " [ string "--left" ; dep lw; string "--right" ; dep rw; string "--pairs_together --PARALLEL_STATS" ]
 
    
let read_normalization seq_type memmory max_cov nb_cpu fastq  : fasta workflow =
	Bistro.Workflow.make ~version:2 ~np:8 ~mem:(99 * 1024) [%sh{|
	TRINITY_PATH=`which Trinity`
	TRINTIY_DIR_PATH=`dirname $TRINITY_PATH`
	READ_NORMALISATION_PATH=$TRINTIY_DIR_PATH/util/insilico_read_normalization.pl
	$READ_NORMALISATION_PATH  {{ config_paired_or_single fastq }} --seqType {{ string seq_type }} --JM {{ seq [ string "$((" ; mem ; string " / 1024))" ]}}G --max_cov {{ int max_cov}} --CPU {{ident np}} --output {{ ident dest }} --no_cleanup
	|}]
	/ selector [ match fastq with
	              | Single_end _ -> "tmp_normalized_reads/single.fa" 
	              | Paired_end _ -> "tmp_normalized_reads/both.fa"
	          ]
	
	
let fastool (fastq : _ fastq workflow) :  fasta workflow =
	Bistro.Workflow.make [%sh {|
	TRINITY_PATH=`which Trinity`
	TRINTIY_DIR_PATH=`dirname $TRINITY_PATH`
	FASTOOL_PATH=$TRINTIY_DIR_PATH/trinity-plugins/fastool/fastool
	$FASTOOL_PATH --illumina-trinity --to-fasta  {{ dep fastq }} > {{ ident dest }}
	|} ]


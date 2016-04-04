open Core_kernel.Std
open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std


let query = "tmp/test/R_sp_Gene_Families/Mus_musculus_F01297.fa"
let fa = "tmp/norm_fasta/Cavia_porcellus.norm.fa"
let db = "test_apytram/db2"
let db_type = "single"

let test_apytram = Apytram.apytram ~plot:true ~i:20 ~tmp:"tmp_apy3" ~fasta:fa ~query:query db_type db

let target_apytram = let open Bistro_app in
                            [[ "test_apytram";"apytram"] %> test_apytram; ]
let _ = Bistro_app.local target_apytram




let seq_dispatcher ?tab_by_family_dir query query_species target seq2sp : fasta workflow = 
       workflow [
            cmd "../bin/SeqDispatcher.py"  [ 
              option (opt "-tab_out_by_family" string) tab_by_family_dir;
			  opt "-q" string query ;
			  opt "-qs" string query_species ;
			  opt "-t" string target ;
			  opt "-t2f" string seq2sp;
			  opt "-out" seq [ dest ; string ("/Trinity." ^ query_species )] ;
            ]
    ]


let query = "output/example_2/Trinity_assembly/Trinity_assembly.Cavia_porcellus_big.fa" (* trinity assembly *)
let target = "tmp/test/R_sp_transcriptome/Mus_musculus_transcriptome.fa" (* annotated transcriptome *)
let seq2sp = "tmp/test/R_sp_Seq_Fam_link//Mus_musculus_Fam_Seq.fa" (* seq2species files *)
let species = "Cavia_porcellus"
let tab_by_family_dir = "tmp_amalgam/Validated_Sequences2Species/"

let seq_dispatcher_out_dir = seq_dispatcher ~tab_by_family_dir:tab_by_family_dir query species target seq2sp 

let target_seq_dispatcher = let open Bistro_app in
                            [[ "test_seq_dispatcher";"seq_dispatcher"] %> seq_dispatcher_out_dir; ]
let _ = Bistro_app.local target_seq_dispatcher



let seq_integrator
      ?tmp
      ?realign_ali
      ?log
      ali
      fasta_to_add
      seq2sp
	  : fasta workflow = 
       workflow [
            cmd "../bin/SeqIntegrator.py"  [ 
              option (opt "-tmp" string) tmp;
			  opt "-ali" string ali ;
			  opt "-fst" string fasta_to_add ;
			  opt "-s2t" string seq2sp;
			  opt "-out" seq [ dest ; string "/seq_integrator"] ;
            ]
    ]

let ali = "example_2/Alignment_data/F01297.fasta"
let fasta_to_add = "test_seq_dispatcher/seq_dispatcher/Trinity."^ species ^"_F01297.fasta"
let seq2sp = "tmp_amalgam/Validated_Sequences2Species/F01297.seq2sp.txt"

let seq_integrator_out_dir = seq_integrator ali fasta_to_add seq2sp

let target_seq_integrator = let open Bistro_app in
                            [[ "test_seq_integrator";"seq_integrator"] %> seq_integrator_out_dir; ]
let _ = Bistro_app.local target_seq_integrator


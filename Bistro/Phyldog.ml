open Core.Std
open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std

type phyldog_configuration = [`phyldog_configuration] directory

type phylotree

let phyldog
    ?datatype
	?dataformat
	?treefile
	?topospecies
	?dlopt
	?equgenomes
	?topogene
	?timelimit
	~threads
	~linkdir
	(seqdir: fasta directory workflow)
	: phylotree directory workflow =
	
	let config_dir = dest // "Configuration" in
	let results_species = dest // "Species_tree/" in
	let results_genes = dest // "Gene_trees/" in
	workflow ~version:4 ~np:threads [
	mkdir_p config_dir;
	mkdir_p results_species;
	mkdir_p results_genes;
	(* Preparing phyldog configuration files*)
	cmd "../bin/PhyldogPrepData.py" [
              option (opt "-datatype" string) datatype ;
              option (opt "-dataformat" string) dataformat ;
              option (opt "-species_tree_file" string) treefile ;
              option (flag string "-topospecies") topospecies ;
              option (opt "-dlopt" string) dlopt ;
              option (opt "-timelimit" int) timelimit ;            
              option (flag string "-equgenomes") equgenomes ;
              option (flag string "-topogene") topogene ;
              opt "-linkdir" dep linkdir;
              opt "-seqdir" dep seqdir;
              opt "-species_tree_resdir" ident results_species;
              opt "-gene_trees_resdir" ident results_genes; 
	  		  opt "-optdir" seq [ ident config_dir ] ;
	  		  ];
	 (* Run phyldog *)
    cmd "mpirun" [
	  	    opt "-np" int threads ;
	  	    string "phyldog";
	  	    seq ~sep:"=" [string "param";  ident (config_dir // "GeneralOptions.txt") ];
	        ];
	]

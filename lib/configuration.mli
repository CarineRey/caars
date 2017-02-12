open Commons

type rna_sample = {
  id : string ;
  species : string ;
  ref_species : string ;
  sample_fastq : string sample_fastq ;
  run_trinity : bool ;
  run_transdecoder : bool ;
  path_assembly : string ;
  given_assembly : bool ;
  run_apytram : bool ;
}

type config_rna_seq = rna_sample list

type t = {
  config_rna_seq : config_rna_seq ;
  apytram_samples: rna_sample list ;
  trinity_samples : rna_sample list ;
  all_ref_samples : rna_sample list ;
  all_ref_species : string list ;
  all_apytram_ref_species : string list ;
  families : string list;
  sample_sheet : string ;
  species_tree_file : string ;
  alignments_dir : string ;
  seq2sp_dir : string ;
  outdir : string ;
  threads : int;
  memory : int;
  run_reconciliation : bool;
}

val load : 
  sample_sheet:string ->
  species_tree_file:string ->
  alignments_dir:string ->
  seq2sp_dir:string ->
  np:int ->
  memory:int ->
  run_reconciliation:bool ->
  outdir:string ->
  t

open Bistro
open Defs
open Wutils

type t = {
  run_reconciliation : bool ;
  configuration_directory : Configuration_directory.t ;
  checked_used_families_all_together : text file ;
  fasta_reads : fasta file OSE_or_PE.t Rna_sample.assoc ;
  normalized_fasta_reads : fasta file OSE_or_PE.t Rna_sample.assoc ;
  trinity_assemblies : fasta file Rna_sample.assoc ;
  trinity_orfs : fasta file Rna_sample.assoc ;
  trinity_assemblies_stats : text file Rna_sample.assoc ;
  trinity_orfs_stats : text file Rna_sample.assoc ;
  trinity_annotated_fams : [`seq_dispatcher] directory Rna_sample.assoc ;
  ref_blast_dbs : blast_db file assoc ;
  reads_blast_dbs : Apytram.compressed_read_db Rna_sample.assoc ;
  apytram_orfs_ref_fams : (Family.t * (Rna_sample.t * Family.t * fasta file) list) list ;
  apytram_checked_families : (Family.t * (Rna_sample.t * Family.t * fasta file) list) list ;
  apytram_annotated_families : (Family.t * fasta file) list ;
  merged_families : (Family.t * [ `seq_integrator ] directory * [ `seq_integrator ] directory option) list ;
  merged_and_reconciled_families : (Family.t * Generax.phylotree directory * [ `seq_integrator ] directory) list ;
  merged_reconciled_and_realigned_families_dirs : [`merged_families_distributor] directory ;
  reconstructed_sequences : [`reconstructed_sequences] directory option ;
  orthologs_per_seq : [`extract_orthologs] directory ;
  final_plots : [`final_plots] directory ;
}

val make :
  ?memory:int ->
  ?nthreads:int ->
  merge_criterion:merge_criterion ->
  filter_threshold:float ->
  refine_ali:bool ->
  run_reconciliation:bool ->
  Dataset.t ->
  t

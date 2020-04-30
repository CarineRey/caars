open Core_kernel
open Bistro
open Bistro_utils.Repo

let ( ++ ) t h = h :: t
let ( +? ) t maybe_h =
  match maybe_h with
  | None -> t
  | Some h -> h :: t

let just_parse_repo (p : Pipeline.t) =
  [
    item ["FamilyMetadata.txt"] (Configuration_directory.family_metadata p.configuration_directory) ;
    item ["SpeciesMetadata.txt"] (Configuration_directory.species_metadata p.configuration_directory) ;
    item ["UsableFamilies.txt"] (Configuration_directory.usable_families p.configuration_directory) ;
    item ["DetectedFamilies.txt"] (Configuration_directory.detected_families p.configuration_directory) ;
    item ["UsedFamilies.txt"] p.checked_used_families_all_together
  ]

let trinity_assemblies (p : Pipeline.t) =
  List.filter_map p.trinity_assemblies ~f:(fun (s, trinity_assembly) ->
      match s.precomputed_assembly with
      | Some _ -> None
      | None ->
        item ["draft_assemblies" ; "raw_assemblies" ; "Draft_assemblies." ^ s.id ^ "_" ^ s.species ^ ".fa"] trinity_assembly
        |> Option.some
    )

let trinity_orfs (p : Pipeline.t) =
  List.filter_map p.trinity_orfs ~f:(fun (s, trinity_orf) ->
      match s.precomputed_assembly with
      | Some _ -> None
      | None ->
        item ["draft_assemblies" ; "cds" ; "Draft_assemblies.cds." ^ s.id ^ "_" ^ s.species ^ ".fa"] trinity_orf
        |> Option.some
    )

let target_to_sample_fasta (s : Rna_sample.t) d = function
  | OSE_or_PE.Single_end se -> [ item [ d ; s.id ^ "_" ^ s.species ^ ".fa" ] se.reads ]
  | Paired_end pe -> [
      item [ d ; s.id ^ "_" ^ s.species ^ ".left.fa" ] pe.reads1 ;
      item [ d ; s.id ^ "_" ^ s.species ^ ".right.fa" ] pe.reads2 ;
    ]

let sample_fasta_repo (p : Pipeline.t) =
  List.concat [
    List.concat_map p.fasta_reads ~f:(fun (s, sample_fasta) -> target_to_sample_fasta s "rna_seq/raw_fasta" sample_fasta) ;
    List.concat_map p.normalized_fasta_reads ~f:(fun (s,norm_fasta) -> target_to_sample_fasta s "rna_seq/norm_fasta" norm_fasta) ;
  ]

let debug_repo (p : Pipeline.t) =
  List.concat [
    List.map p.trinity_assemblies_stats ~f:(fun (s, trinity_assembly_stats) ->
        item ["debug" ; "trinity_assembly" ; "trinity_assemblies_stats" ; "Trinity_assemblies." ^ s.id ^ "_" ^ s.species ^ ".stats"] trinity_assembly_stats
      ) ;
    List.map p.trinity_orfs_stats ~f:(fun (s, trinity_orfs_stats) ->
        item ["debug" ; "trinity_assembly" ; "trinity_assemblies_stats" ; "Transdecoder_cds." ^ s.id ^ "_" ^ s.species ^ ".stats"] trinity_orfs_stats
      ) ;
    List.map p.trinity_annotated_fams ~f:(fun (s, trinity_annotated_fams) ->
        item ["debug" ; "trinity_blast_annotation" ; "trinity_annotated_fams" ; s.id ^ "_" ^ s.species ^ ".vs." ^ String.concat ~sep:"_" s.reference_species] trinity_annotated_fams
      ) ;
    List.map p.ref_blast_dbs ~f:(fun (ref_species, blast_db) ->
        [ "debug" ; "trinity_blast_annotation" ; "ref_blast_db" ; ref_species ] %> blast_db
      ) ;
    List.map p.reads_blast_dbs ~f:(fun (s,blast_db) ->
        [ "debug" ; "rna_seq" ;"rep_cluster_blast_db" ; s.id ^ "_" ^ s.species ] %> blast_db.cluster_rep_blast_db
      )
    ;
    List.map p.apytram_orfs_ref_fams ~f:(fun (_, fws) ->
        List.map fws ~f:(fun (s, fam, apytram_result) ->
            [ "debug" ; "apytram_assembly" ; "apytram_results_by_ref_by_group_by_fam" ; fam.name ; s.id ^ "_" ^ s.species ^ ".fa" ] %> apytram_result
          )
      ) |> List.concat
    ;
    List.map p.apytram_checked_families ~f:(fun (_, fws) ->
        List.map fws ~f:(fun (s, fam, apytram_result) ->
            [ "debug" ; "apytram_assembly" ; "apytram_checked_families" ; fam.name ; s.id ^ "_" ^ s.species ^ ".fa"] %> apytram_result
          )
      ) |> List.concat
    ;
    List.map p.apytram_annotated_families ~f:(fun (fam, fw) ->
        [["debug" ; "apytram_assembly" ;"apytram_annotated_sequences"; fam.name ] %> fw]
      ) |> List.concat
    ;
    List.concat_map p.merged_families ~f:(fun (fam, merged_family, merged_and_filtered_family) ->
        match (merged_family, merged_and_filtered_family) with
        | (w1, Some w2) ->  [ [ "debug" ; "merged_families" ; fam.name  ] %> w1; [ "debug" ; "merged_filtered_families" ; fam.name  ] %> w2 ]
        | (w1, None) -> [[ "debug" ; "merged_families" ; fam.name  ] %> w1]
      )
  ]

let full_repo ?(get_reads = false) ?(debug = false) (p : Pipeline.t) =
  List.concat [
    just_parse_repo p ;
    trinity_assemblies p ;
    trinity_orfs p ;
    ([]
     ++ item ["assembly_results_by_fam" ] (Workflow.select p.merged_reconciled_and_realigned_families_dirs ["out/"])
     ++ item ["all_fam.seq2sp.tsv"]       (Workflow.select p.orthologs_per_seq ["all_fam.seq2sp.tsv"])
     ++ item ["plots"] p.final_plots
     +? Option.map p.reconstructed_sequences ~f:(fun w -> item ["assembly_results_only_seq"] (Workflow.select w ["assemblies/"]))
     +? (
       if p.run_reconciliation then
         Some (item ["all_fam.orthologs.tsv"] (Workflow.select p.orthologs_per_seq ["all_fam.orthologs.tsv"]))
       else None
     )
    ) ;
    if get_reads then sample_fasta_repo p else [] ;
    if debug then debug_repo p else [] ;
  ]

let precious_repo (p : Pipeline.t) =
  let pi = precious_item in
  let normalized_fasta_samples =
    List.concat_map p.normalized_fasta_reads ~f:(function
        | (_, Single_end se) -> [ pi se.reads ]
        | (_, Paired_end pe) -> [ pi pe.reads1 ; pi pe.reads2 ]
      )
  in
  let reads_blast_dbs =
    List.concat_map p.reads_blast_dbs ~f:(fun (_, x) ->
        [pi x.concat_fasta; pi x.index_concat_fasta; pi x.rep_cluster_fasta; pi x.reformated_cluster; pi x.index_cluster ; pi x.cluster_rep_blast_db ]
      )
  in
  let get_merged_families = function
    | (_, w1, Some w2) -> [pi w1; pi w2]
    | (_, w1, None) -> [pi w1]
  in
  List.concat [
    [pi p.configuration_directory];
    normalized_fasta_samples ;
    List.map p.trinity_assemblies ~f:(fun (_, x) -> pi x) ;
    List.map p.trinity_orfs ~f:(fun (_, x) -> pi x) ;
    reads_blast_dbs ;
    List.map p.trinity_annotated_fams ~f:(fun (_, x) -> pi x) ;
    List.concat_map p.merged_families ~f:get_merged_families;
    List.map p.merged_and_reconciled_families ~f:(fun (_, x, _) -> pi x) ;
    List.concat_map p.merged_and_reconciled_families ~f:(fun (_ , w1, w2) -> [pi w1; pi w2;]) ;
    List.concat_map p.apytram_checked_families ~f:(fun (_, fws) -> List.map fws ~f:(fun (_, _, x) -> pi x)) ;
    List.map p.apytram_annotated_families ~f:(fun (_, x) -> pi x) ;
  ]

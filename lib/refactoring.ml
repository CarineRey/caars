open Core_kernel
open Bistro

type 'a assoc = (string, 'a) List.Assoc.t

let ( $ ) a k = List.Assoc.find_exn ~equal:String.equal a k

let assoc keys ~f =
  List.map keys ~f:(fun k -> k, f k)

let assoc_map a ~f =
  List.map a ~f:(fun (k, v) -> k, f k v)

let assoc_filter_map a ~f =
  List.filter_map a ~f:(fun (k, v) -> Option.map (f k v) ~f:(fun v -> k, v))

let assoc_opt keys ~f =
  List.filter_map keys ~f:(fun k ->
      Option.map (f k) ~f:(fun v -> k, v)
    )

let id_concat xs =
  List.filter xs ~f:(String.( <> ) "")
  |> String.concat ~sep:"_"

let descr ?label d =
  match label with
  | None -> d
  | Some label -> sprintf "%s:%s" d label

class type blast_db = object
  inherit text
  method format : [`blast_db]
end

module SE_or_PE_file = struct
  type single_end_orientation =
    | F
    | R
    | US

  type paired_end_orientation =
    | FR
    | RF
    | UP

  type 'a t =
    | Single_end of {
        reads : 'a file ;
        orientation : single_end_orientation ;
      }
    | Paired_end of {
        reads1 : 'a file ;
        reads2 : 'a file ;
        orientation : paired_end_orientation ;
      }

  let se reads orientation = Single_end { reads ; orientation }
  let pe reads1 reads2 orientation = Paired_end { reads1 ; reads2 ; orientation }
end

module Sample_file = struct
  type t =
    | Fastq_sample_file of fastq SE_or_PE_file.t
    | Fasta_sample_file of fasta SE_or_PE_file.t
end

module Trinity = struct
  open Bistro.Shell_dsl
  open Commons
  include Trinity

  let single_stranded_or_unstranded = function
    | SE_or_PE_file.F -> string "--SS_lib_type F"
    | R -> string "--SS_lib_type R"
    | US -> string ""

  let paired_stranded_or_unstranded = function
    | SE_or_PE_file.RF -> string "--SS_lib_type RF"
    | FR -> string "--SS_lib_type FR"
    | UP -> string ""

  let config_trinity_fasta_paired_or_single = function
    | SE_or_PE_file.Single_end se ->
      seq ~sep: " " [ string "--single" ; dep se.reads ; single_stranded_or_unstranded se.orientation ]
    | Paired_end pe ->
      seq ~sep: " " [ string "--left" ; dep pe.reads1 ; string "--right" ; dep pe.reads2 ; paired_stranded_or_unstranded pe.orientation ]

  let fasta_read_normalization_get_output ~fasta ~dest=
    let (vars, code) = match fasta with
      | SE_or_PE_file.Single_end _ -> (["DEST", dest;
                                "SINGLELINK", string "`readlink single.norm.fa`"],
                               {| mv $SINGLELINK $DEST/"single.norm.fa"|})
      | Paired_end _ -> (["DEST", dest;
                                "LEFTLINK", string "`readlink left.norm.fa`";
                                "RIGHTLINK", string "`readlink right.norm.fa`"],
                               {|echo $LEFTLINK ; mv $LEFTLINK $DEST/"left.norm.fa"; mv $RIGHTLINK $DEST/"right.norm.fa"|})
    in
    bash_script vars code

  let fasta_read_normalization
    ?(descr = "")
    max_cov
    ~threads
    ?(memory = 1)
    ?(max_memory = 1)
    (fasta : fasta SE_or_PE_file.t)
    : fasta SE_or_PE_file.t =
  let bistro_memory =
    if max_memory > 2 then Int.(min max_memory (memory * 2))
    else 1
  in
  let given_mem =
    if bistro_memory > 2 then Int.(bistro_memory / 2)
    else 1
  in
  (* reserve more memory by bistro than given to normalization tools*)
  let output_dir =
    Workflow.shell ~descr:("fasta_read_normalization" ^ descr) ~version:2 ~np:threads ~mem:(Workflow.int (1024 * bistro_memory)) [
      mkdir_p dest;
      mkdir_p tmp ;
      within_container img (
        and_list [
          cmd "Trinity" [
            string "--no_version_check";
            opt "--max_memory" ident (seq [ string "$((" ; int given_mem ; string " / 1024))G" ]) ;
            opt "--CPU" ident np ;
            string "--just_normalize_reads";
            opt "--normalize_max_read_cov" int max_cov ;
            config_trinity_fasta_paired_or_single fasta ;
            string "--seqType fa" ;
            opt "--output" seq [ ident tmp ; string "/trinity"] ;
          ];
          cd (tmp // "trinity/insilico_read_normalization") ;
          cmd "sh" [ file_dump (fasta_read_normalization_get_output ~fasta ~dest) ];
        ]
      )
    ]
  in
  match fasta with
  | Single_end se -> SE_or_PE_file.se (Workflow.select output_dir ["single.norm.fa"]) se.orientation
  | Paired_end pe ->
    SE_or_PE_file.pe
      (Workflow.select output_dir ["left.norm.fa"])
      (Workflow.select output_dir ["right.norm.fa"])
      pe.orientation

  let trinity_fasta
      ?label
      ?full_cleanup
      ?no_normalization
      ~threads
      ?(memory = 1)
      (sample_fasta : fasta SE_or_PE_file.t)
    : fasta file =
    Workflow.shell ~descr:(descr ?label "Trinity") ~np:threads ~mem:(Workflow.int (1024 * memory)) [
      mkdir_p dest;
      cmd "Trinity" ~img [
        string "--no_version_check";
        opt "--max_memory" ident (seq [ string "$((" ; mem ; string " / 1024))G" ]) ;
        opt "--CPU" ident np ;
        option (flag string "--full_cleanup") full_cleanup ;
        option (flag string "--no_normalize_reads") no_normalization ;
        config_trinity_fasta_paired_or_single sample_fasta;
        string "--seqType fa" ;
        opt "--output" seq [ ident dest ; string "/trinity"] ;
      ];
      cmd "sed" [
        string "-re";
        string {|"s/(>[_a-zA-Z0-9]*)( len=[0-9]* path=.*)/\1/"|};
        string "-i";
        seq [ident dest; string "/trinity.Trinity.fasta";];
      ];
    ]
    |> Fn.flip Workflow.select [ "trinity.Trinity.fasta" ]
end

module Rna_sample = struct
  type t = {
    id : string ;
    group_id : string ;
    species : string ;
    reference_species : string list ;
    sample_file : Sample_file.t ;
    run_trinity : bool ;
    run_transdecoder : bool ;
    run_apytram : bool ;
    precomputed_assembly : string option ;
  }
  type 'a assoc = (t * 'a) list
end

module Configuration = struct
  type merge_criterion =
    | Merge
    | Length
    | Length_complete

  type t = {
    samples : Rna_sample.t list ;
    reference_species : string list ;
    reference_transcriptome : fasta file assoc ;
    memory : int ;
    nthreads : int ;
  }

  let reference_transcriptome config spe = config.reference_transcriptome $ spe

  let samples_with_trinity_run config =
    List.filter config.samples ~f:(fun s -> s.run_trinity)

  let make ?(nthreads = 1) ?(memory = 1) ~reference_transcriptome ~samples = {
    samples ;
    reference_species = (
      List.concat_map samples ~f:(fun (it : Rna_sample.t) -> it.reference_species)
      |> List.dedup_and_sort ~compare:String.compare
    ) ;
    reference_transcriptome ;
    nthreads ; memory ;
  }
end

module Pipeline = struct
  type t = {
    fasta_reads : fasta SE_or_PE_file.t Rna_sample.assoc ;
    normalized_fasta_reads : fasta SE_or_PE_file.t Rna_sample.assoc ;
    trinity_assemblies : fasta file Rna_sample.assoc ;
    trinity_orfs : fasta file Rna_sample.assoc ;
    trinity_assemblies_stats : text file Rna_sample.assoc ;
    trinity_orfs_stats : text file Rna_sample.assoc ;
    ref_blast_dbs : blast_db file assoc ;
  }

  let rna_sample_needs_rna (s : Rna_sample.t) =
    match s.run_apytram, s.run_trinity, s.precomputed_assembly with
    | true, _, _          -> true
    | false, true, Some _ -> false
    | false, true, None   -> true
    | false, false, _     -> false

  let fasta_reads (config : Configuration.t) =
    assoc_opt config.samples ~f:(fun s ->
        if rna_sample_needs_rna s then
          let open SE_or_PE_file in
          let fq2fa ?(label = "") x =
            let descr = id_concat [s.id ; s.species ; label] in
            Trinity.fastq2fasta ~descr x
          in
          match s.sample_file with
          | Fasta_sample_file fa -> Some fa
          | Fastq_sample_file (Single_end se) ->
            Some (SE_or_PE_file.se (fq2fa se.reads) se.orientation)
          | Fastq_sample_file (Paired_end pe) ->
            Some (SE_or_PE_file.pe
                    (fq2fa ~label:"left" pe.reads1)
                    (fq2fa ~label:"right" pe.reads2)
                    pe.orientation)
        else None
      )

  let normalize_fasta_reads fasta_reads memory max_memory threads =
     assoc_map fasta_reads ~f:(fun (s : Rna_sample.t) fa ->
        let max_cov = 20 in
        let descr = id_concat [s.id ; s.species] in
        Trinity.fasta_read_normalization ~descr max_cov ~threads ~memory ~max_memory fa
      )

  let trinity_assemblies_of_norm_fasta normalized_fasta_reads ~memory ~nthreads =
    assoc_filter_map normalized_fasta_reads ~f:(fun (s : Rna_sample.t) normalized_fasta_reads ->
        match s.run_trinity, s.precomputed_assembly with
        | true, None ->
          let label = id_concat [s.id ; s.species] in
          Trinity.trinity_fasta ~label ~no_normalization:true ~full_cleanup:true ~memory ~threads:nthreads normalized_fasta_reads
          |> Option.some
        | _, Some assembly_path -> Some (Workflow.input assembly_path)
        | (_, _)   -> None
      )

  let transdecoder_orfs_of_trinity_assemblies trinity_assemblies ~memory ~nthreads =
    assoc_map trinity_assemblies ~f:(fun (s : Rna_sample.t) trinity_assembly ->
        match s.run_transdecoder, s.precomputed_assembly with
        | true, None ->
          let pep_min_length = 50 in
          let retain_long_orfs = 150 in
          let descr = id_concat ["Assembly" ; s.id ; s.species] in
          Transdecoder.transdecoder
            ~descr ~retain_long_orfs ~pep_min_length ~only_best_orf:false
            ~memory ~threads:nthreads
            trinity_assembly
        | (false, _ )
        | (true, Some _) -> trinity_assembly
      )

  let assemblies_stats_of_assemblies assemblies =
    assoc_filter_map assemblies ~f:(fun (s : Rna_sample.t) assembly ->
        match s.precomputed_assembly with
        | Some _ -> None
        | None ->
          Some (Trinity.assembly_stats ~descr:(s.id ^ "_" ^ s.species) assembly)
      )

  let ref_blast_dbs (config : Configuration.t) =
    assoc config.reference_species ~f:(fun ref_species ->
        let fasta = Configuration.reference_transcriptome config ref_species in
        let parse_seqids = true in
        let dbtype = "nucl" in
        BlastPlus.makeblastdb ~parse_seqids ~dbtype  ("DB_" ^ ref_species) fasta
      )

  let make (config : Configuration.t) =
    let memory_per_sample, threads_per_sample =
      let nb_samples = List.length config.samples in
      Int.(max 1 (config.memory / (max 1 nb_samples))), Stdlib.(max 1 (config.nthreads / (max 1 nb_samples)))
    in
    (* let memory_per_thread = Int.(max 1 (config.memory / config.nthreads)) in *)
    let ref_blast_dbs = ref_blast_dbs config in
    let fasta_reads = fasta_reads config in
    let normalized_fasta_reads = normalize_fasta_reads fasta_reads memory_per_sample config.memory threads_per_sample in
    let trinity_assemblies = trinity_assemblies_of_norm_fasta normalized_fasta_reads ~memory:memory_per_sample ~nthreads:threads_per_sample in
    let trinity_orfs = transdecoder_orfs_of_trinity_assemblies trinity_assemblies ~memory:memory_per_sample ~nthreads:threads_per_sample in
    let trinity_assemblies_stats = assemblies_stats_of_assemblies trinity_assemblies in
    let trinity_orfs_stats = assemblies_stats_of_assemblies trinity_orfs in
    { ref_blast_dbs ; fasta_reads ; normalized_fasta_reads ;
      trinity_assemblies ; trinity_orfs ; trinity_assemblies_stats ;
      trinity_orfs_stats }
end

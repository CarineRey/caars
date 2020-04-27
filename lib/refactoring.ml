open Core_kernel
open Bistro

type 'a assoc = (string, 'a) List.Assoc.t

let ( $ ) a k = List.Assoc.find_exn ~equal:String.equal a k

let assoc_opt keys ~f =
  List.filter_map keys ~f:(fun k ->
      Option.map (f k) ~f:(fun v -> k, v)
    )

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
  }

  let reference_transcriptome config spe = config.reference_transcriptome $ spe

  let make ~reference_transcriptome ~samples = {
    samples ;
    reference_species = (
      List.concat_map samples ~f:(fun (it : Rna_sample.t) -> it.reference_species)
      |> List.dedup_and_sort ~compare:String.compare
    ) ;
    reference_transcriptome ;
  }
end

module Pipeline = struct
  type t = {
    fasta_reads : (Rna_sample.t, fasta SE_or_PE_file.t) List.Assoc.t ;
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
          let fq2fa ?label x =
            let label = Option.value_map label ~default:"" ~f:(( ^ ) "_") in
            let descr = sprintf "%s_%s%s" s.id s.species label in
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

  let ref_blast_dbs (config : Configuration.t) =
    List.map config.reference_species ~f:(fun ref_species ->
        let fasta = Configuration.reference_transcriptome config ref_species in
        let parse_seqids = true in
        let dbtype = "nucl" in
        (ref_species, BlastPlus.makeblastdb ~parse_seqids ~dbtype  ("DB_" ^ ref_species) fasta)
      )

  let make configuration =
    let ref_blast_dbs = ref_blast_dbs configuration in
    let fasta_reads = fasta_reads configuration in
    { ref_blast_dbs ; fasta_reads }
end

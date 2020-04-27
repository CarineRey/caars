open Core_kernel
open Bistro

type 'a assoc = (string, 'a) List.Assoc.t

let ( $ ) a k = List.Assoc.find_exn ~equal:String.equal a k

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
end

module Sample_file = struct
  type t =
    | Fastq_sample_file of fastq SE_or_PE_file.t
    | Fasta_sample_file of fasta SE_or_PE_file.t
end

module Rna_sample = struct
  type t = {
    id : string ;
    species : string ;
    apytram_group : string ;
    ref_species : string list;
    sample_file : Sample_file.t ;
    run_trinity : bool ;
    run_transdecoder : bool ;
    path_assembly : string ;
    given_assembly : bool ;
    run_apytram : bool ;
  }
end

module Configuration = struct
  type merge_criterion =
    | Merge
    | Length
    | Length_complete

  type item = {
    id : string ;
    group_id : string ;
    species : string ;
    reference_species : string ;
    sample : Rna_sample.t ;
  }

  type t = {
    samples : item list ;
    reference_species : string list ;
    reference_transcriptome : fasta file assoc ;
  }

  let reference_transcriptome config spe = config.reference_transcriptome $ spe

  let make ~reference_transcriptome ~samples = {
    samples ;
    reference_species = (
      List.map samples ~f:(fun (it : item) -> it.reference_species)
      |> List.dedup_and_sort ~compare:String.compare
    ) ;
    reference_transcriptome ;
  }
end

module Pipeline = struct
  type t = {
    ref_blast_dbs : blast_db file assoc ;
  }

  let ref_blast_dbs (config : Configuration.t) =
    List.map config.reference_species ~f:(fun ref_species ->
        let fasta = Configuration.reference_transcriptome config ref_species in
        let parse_seqids = true in
        let dbtype = "nucl" in
        (ref_species, BlastPlus.makeblastdb ~parse_seqids ~dbtype  ("DB_" ^ ref_species) fasta)
      )

  let make configuration =
    let ref_blast_dbs = ref_blast_dbs configuration in
    { ref_blast_dbs }
end

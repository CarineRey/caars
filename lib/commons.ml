(*
# File: Commons.ml
# Created by: Carine Rey
# Created on: March 2016
#
#
# Copyright 2016 Carine Rey
# This software is a computer program whose purpose is to assembly
# sequences from RNA-Seq data (paired-end or single-end) using one or
# more reference homologous sequences.
# This software is governed by the CeCILL license under French law and
# abiding by the rules of distribution of free software.  You can  use,
# modify and/ or redistribute the software under the terms of the CeCILL
# license as circulated by CEA, CNRS and INRIA at the following URL
# "http://www.cecill.info".
# As a counterpart to the access to the source code and  rights to copy,
# modify and redistribute granted by the license, users are provided only
# with a limited warranty  and the software's author,  the holder of the
# economic rights,  and the successive licensors  have only  limited
# liability.
# In this respect, the user's attention is drawn to the risks associated
# with loading,  using,  modifying and/or developing or reproducing the
# software by the user in light of its specific status of free software,
# that may mean  that it is complicated to manipulate,  and  that  also
# therefore means  that it is reserved for developers  and  experienced
# professionals having in-depth computer knowledge. Users are therefore
# encouraged to load and test the software's suitability as regards their
# requirements in conditions enabling the security of their systems and/or
# data to be ensured and,  more generally, to use and operate it in the
# same conditions as regards security.
# The fact that you are presently reading this means that you have had
# knowledge of the CeCILL license and that you accept its terms.
*)

open Bistro.Std
open Bistro.EDSL
open Bistro_bioinfo.Std

type ('a,'b) either =
  | Left of 'a
  | Right of 'b

type ('a,'b) fastX =
  | Fasta of 'a
  | Fastq of 'b

type orientation_single =
  | F
  | R
  | US

type orientation_paired =
  | FR
  | RF
  | UP

type ('a,'b) sample_file =
  | Sample_fastq of 'a
  | Sample_fasta of 'b

type 'a sample_fastq =
  | Fastq_Single_end of 'a * orientation_single
  | Fastq_Paired_end of 'a * 'a * orientation_paired

type 'a sample_fasta =
  | Fasta_Single_end of 'a * orientation_single
  | Fasta_Paired_end of 'a * 'a * orientation_paired

let sample_fastq_map f = function
  | Fastq_Single_end ( x , o ) ->  Fastq_Single_end ( f x , o )
  | Fastq_Paired_end ( lx, rx, o ) -> Fastq_Paired_end (f lx, f rx, o)

let sample_fasta_is_paired = function
  | Fasta_Single_end _ ->  false
  | Fasta_Paired_end _ -> true

let sample_fastq_orientation = function
  | Fastq_Single_end ( x , o ) ->  Left o
  | Fastq_Paired_end ( lx, rx, o ) -> Right o

let sample_fasta_map f = function
  | Fasta_Single_end ( x , o ) ->  Fasta_Single_end ( f x , o )
  | Fasta_Paired_end ( lx, rx, o ) -> Fasta_Paired_end (f lx, f rx, o)

let sample_file_map f = function
  | Sample_fasta x ->  Sample_fasta (sample_fasta_map f x)
  | Sample_fastq x ->  Sample_fastq (sample_fastq_map f x)

let sample_fastq_is_paired = function
  | Fastq_Single_end _ ->  false
  | Fastq_Paired_end _ -> true

let sample_fasta_orientation = function
  | Fasta_Single_end ( x , o ) ->  Left o
  | Fasta_Paired_end ( lx, rx, o ) -> Right o


type rna_sample = {
  id : string ;
  species : string ;
  ref_species : string list;
  sample_file : ( string sample_fastq, string sample_fasta ) sample_file ;
  sample_fastq: string sample_fastq ;
  run_trinity : bool ;
  run_transdecoder : bool ;
  path_assembly : string ;
  given_assembly : bool ;
  run_apytram : bool ;
}

type config_rna_seq = rna_sample list
type output = [ `amalgam_output ]
type configuration_dir = [ `configuration ]

type sp2seq_link
type tabular
type index
type cdhit
type blast_db = [`blast_db] directory


type compressed_read_db = {
  s : rna_sample;
  concat_fasta : fasta workflow;
  index_concat_fasta : index workflow;
  rep_cluster_fasta : fasta workflow;
  reformated_cluster : fasta workflow;
  index_cluster : index workflow;
  cluster_rep_blast_db : blast_db workflow;
}




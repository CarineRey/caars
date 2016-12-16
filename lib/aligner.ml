open Core.Std
open Bistro.Std
open Bistro.EDSL
open Bistro_bioinfo.Std


let mafft
  ?treein
  ~auto
  ?maxiterate
  ~threads
  (fa:fasta workflow) : fasta workflow =
  workflow ~descr:"mafft" ~np:threads [
    cmd "mafft" ~stdout:dest [
      flag string "--auto" auto;
      option (opt "--treein" dep) treein;
      option (opt "--maxiterate" int) maxiterate;
      dep fa ;
    ]
  ]

open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std


let mafft 
  ?treein
  ?(auto = true)
  ?maxiterate
  (fa:fasta workflow) : fasta workflow =
  workflow ~descr:"mafft" ~np [
    cmd "mafft" ~stdout:dest [
      option (flag string "--auto") auto;
      option (opt "--treein" dep) treein;
      option (opt "--maxiterate" int) maxiterate;
      dep fa ;
    ]
  ]

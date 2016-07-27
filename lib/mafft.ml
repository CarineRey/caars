open Bistro.Std
open Bistro.EDSL_sh
open Bistro_bioinfo.Std


let mafft (fa:fasta workflow) : fasta workflow =
  workflow ~descr:"mafft" [
    cmd "mafft" ~stdout:dest [
      dep fa ;
    ]
  ]

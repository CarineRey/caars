open Core.Std
open Bistro.Std
open Bistro.EDSL
open Bistro_bioinfo.Std


let mafft
  ?(descr="")
  ~treein
  ~auto
  ?maxiterate
  ~threads
  (fa:fasta workflow) : fasta workflow =
   let script = [%bistro {|
       cp {{ dep treein }} tmp.tree
       grep "^>" {{ dep fa }} | sed "s/>//" | awk '{ print $0 ":" NR }' > tmp.num
       for i in $(cat tmp.num)
       do
         s=`echo $i | cut -f 1 -d ":"`
         n=`echo $i | cut -f 2 -d ":"`
         sed "s/$s/$n/" tmp.tree -i
       done
       nw2nhx.py tmp.tree > tmp.tree
       newick2mafft.rb tmp.tree > tmp.mafft
       mafft --treein tmp.mafft {{ dep fa }} > {{ ident dest }}
       |} ]
  in
  workflow ~descr:("mafft"^descr) ~version:2 ~np:threads [
  cd tmp;
  cmd "sh" [ file_dump script ];
  ]


(*  workflow ~descr:"mafft" ~np:threads [
    cd tmp;
    let w =  match treein with
      | None   -> cmd "mafft" ~stdout:dest [
                      flag string "--auto" auto;
                      option (opt "--treein" dep) treein;
                      option (opt "--maxiterate" int) maxiterate;
                      dep fa ;
                      ]
     | Some _ -> let script = [%bistro {|
       cp {{ dep treein }} tmp.tree
       grep "^>" {{ dep fa }} | sed "s/>//" | awk '{ print $0 ":" NR }' > tmp.num
       for i in $(cat tmp.num)
       do
       s=`echo $i | cut -f 1 -d ":"`
       n=`echo $i | cut -f 2 -d ":"`
       sed "s/$s/$n/" tmp.tree -i
       done
       newick2mafft.rb tmp.tree > tmp.mafft
        mafft --treein tmp.mafft {{ dep fa }} > {{ dep dest }}
       |} ]
       in
       cmd "sh" [ file_dump script ]
    in
    w
   ;
  ]
*)

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
       nw2nhx.py tmp.tree tmp.tree
       newick2mafft.rb tmp.tree > tmp.mafft
       mafft --treein tmp.mafft {{ dep fa }} > {{ ident dest }}
       |} ]
  in
  workflow ~descr:("mafft"^descr) ~version:2 ~np:threads [
  cd tmp;
  cmd "sh" [ file_dump script ];
  ]


let mafft_from_nogap
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
       nw2nhx.py tmp.tree tmp.tree
       newick2mafft.rb tmp.tree > tmp.mafft
       awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' {{ dep fa }} | sed "s/-//g" > tmp.nogap.fa
       diff tmp.nogap.fa  {{ dep fa }}
       mafft --treein tmp.mafft tmp.nogap.fa > {{ident dest }}
       |} ]
  in
  workflow ~descr:("mafft_from_no_gap"^descr) ~version:4 ~np:threads [
  cd tmp;
  cmd "sh" [ file_dump script ];
  ]


let muscle
    ?(descr="")
    ~maxiters fa: fasta workflow =
    workflow ~descr:("muscle"^descr) ~version:1 ~np:1 [
    mkdir_p tmp;
    cmd "muscle" [
        opt "-in" dep fa ;
        opt "-out" ident dest ;
        opt "-maxiters" int maxiters;
        ];
    ]
let muscletreein
    ?(descr="")
    ~treein
    ~maxiters fa: fasta workflow =
    let treenw = tmp // "treenw.nw" in  
    let script_nhx_nw = [%bistro {|
      nhx2nw.py {{dep treein}} {{ident treenw}}|}]
    in
    workflow ~descr:("muscletreein"^descr) ~version:1 ~np:1 [
    mkdir_p tmp;
    cmd "sh" [ file_dump script_nhx_nw ];
    cmd "muscle" [
        opt "-in" dep fa ;
        opt "-out" ident dest ;
        opt "-usetree" ident treenw;
        opt "-maxiters" int maxiters;
        ];
    ]
let musclenogap
    ?(descr="")
    ?(treein)
    ~maxiters fa : fasta workflow =
    let nogapfa = tmp // "nogap.fa" in  
    let scriptnogap = [%bistro {|
       awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' {{ dep fa }} | sed "s/-//g" > {{ident nogapfa}} |} ]
    in
    
    workflow ~descr:("musclenogap"^descr) ~version:1 ~np:1 [
    
    mkdir_p tmp;
    cmd "sh" [ file_dump scriptnogap ];
    cmd "muscle" [
        opt "-in" ident nogapfa ;
        opt "-out" ident dest ;
        opt "-maxiters" int maxiters;
        ];
    ]
let musclenogaptreein
    ?(descr="")
    ~treein
    ~maxiters fa : fasta workflow =
    let nogapfa = tmp // "nogap.fa" in  
    let scriptnogap = [%bistro {|
       awk '!/^>/ { printf "%s", $0; n = "\n" } /^>/ { print n $0; n = "" } END { printf "%s", n } ' {{ dep fa }} | sed "s/-//g" > {{ident nogapfa}} |} ]
    in
    let treenw = tmp // "treenw.nw" in  
    let script_nhx_nw = [%bistro {|
      nhx2nw.py {{dep treein}} {{ident treenw}}|}]
    in
    
    workflow ~descr:("musclenogaptreein"^descr) ~version:1 ~np:1 [
    
    mkdir_p tmp;

    cmd "sh" [ file_dump script_nhx_nw ];
    cmd "sh" [ file_dump scriptnogap ];
    cmd "muscle" [
        opt "-in" ident nogapfa ;
        opt "-out" ident dest ;
        opt "-usetree" ident treenw;
        opt "-maxiters" int maxiters;
        ];
    ]

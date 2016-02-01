ocamlbuild -use-ocamlfind -tag thread -pkgs bistro.utils,bistro.bioinfo Main.byte && ./Main.byte  example/conf.tsv example/ref_seqs.tsv

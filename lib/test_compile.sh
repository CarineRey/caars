ocamlbuild -use-ocamlfind -tag thread -pkgs bistro.utils,bistro.bioinfo,ppx_bistro Test_bistro.byte && ./Test_bistro.byte test.fasta test_out

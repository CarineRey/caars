amalgam:
	ocamlbuild -use-ocamlfind -I lib -tag thread -pkgs bistro.utils,bistro.bioinfo,ppx_bistro app/amalgam.byte

test:
	ocamlbuild -use-ocamlfind -I lib -tag thread -pkgs bistro.utils,bistro.bioinfo,ppx_bistro app/test_bistro.byte
	app/test_bistro.byte data/test.fasta test_out

clean:
	ocamlbuild -clean
	rm -f utils/lib/*.pyc

.PHONY: amalgam test clean

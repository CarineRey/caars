amalgam:
	ocamlbuild -use-ocamlfind -I lib -tag thread -pkgs bistro.utils,bistro.bioinfo,ppx_bistro app/amalgam.byte

.PHONY: amalgam clean

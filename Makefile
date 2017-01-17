amalgam:
	ocamlbuild -use-ocamlfind -I lib -tag thread -pkgs bistro.utils,bistro.bioinfo,ppx_bistro app/amalgam_app.byte

test:
	cd example && bash Launch_amalgam.sh

clean_test:
	cd example && rm -r working_dir/ output_dir/

clean:
	ocamlbuild -clean
	rm -f utils/lib/*.pyc

.PHONY: amalgam test clean_test clean

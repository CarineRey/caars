caars:
	ocamlbuild -use-ocamlfind -I lib -tag thread -pkgs bistro.utils,bistro.bioinfo app/caars_app.byte
	mv caars_app.byte caars
test:
	cd example && bash Launch_CAARS.sh

test2:
	cd example && bash Launch_CAARS2.sh

clean_test:
	cd example && rm -r working_dir/ output_dir/

clean_test2:
	cd example && rm -r working2_dir/ output2_dir/

clean:
	ocamlbuild -clean
	rm -f utils/lib/*.pyc

.PHONY: caars test clean_test clean

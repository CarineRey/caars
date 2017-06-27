caars:
	jbuilder build
	cp _build/default/app/caars_app.exe caars
test:
	cd example && bash Launch_CAARS.sh

test2:
	cd example && bash Launch_CAARS2.sh

clean_test:
	cd example && (rm -r working_dir/ output_dir/ || exit 0)

clean_test2:
	cd example && (rm -r working2_dir/ output2_dir/ || exit 0)

clean:
	ocamlbuild -clean
	rm -f utils/lib/*.pyc

.PHONY: caars test clean_test clean

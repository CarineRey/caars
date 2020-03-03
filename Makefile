caars:
	dune build
	cp _build/default/app/caars_app.exe caars
test:
	cd example && bash Launch_CAARS.sh
test_docker:
	cd example && bash Launch_CAARS.sh docker

test2:
	cd example && bash Launch_CAARS2.sh

test2_docker:
	cd example && bash Launch_CAARS2.sh docker

clean_test:
	cd example && (rm -r working_dir/ output_dir/ || exit 0)

clean_test2:
	cd example && (rm -r working2_dir/ output2_dir/ || exit 0)

clean:
	ocamlbuild -clean
	rm -f utils/lib/*.pyc

build_caars_env_docker:
	cd etc && ./build_caars_env.sh
build_caars_docker:
	cd etc && ./build_caars.sh
build_caars_dev_docker:
	cd etc && ./build_caars_dev.sh
push_caars_env_docker:
	cd etc && ./build_caars_env.sh push_yes
push_caars_docker:
	cd etc && ./build_caars.sh push_yes
push_caars_dev_docker:
	cd etc && ./build_caars_dev.sh push_yes

.PHONY: caars test clean_test clean test2 clean_test2 build_caars_env_docker build_caars_docker build_caars_dev_docker

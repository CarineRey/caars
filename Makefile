##########################
# build caars executable #
##########################

caars:
	dune build
	cp _build/default/app/caars_app.exe caars

clean:
	dune clean
	rm -f utils/lib/*.pyc


################################
# build and push caars dockers #
################################

# caars_env is a docker containing all caars dependencies
# caars is a docker containing all caars dependencies and caars

build_caars_env_docker:
	cd etc && ./build_caars_env.sh
build_caars_docker:
	cd etc && ./build_caars_docker.sh
build_caars_docker_travis:
	cd etc && ./build_caars_docker.sh build ${TRAVIS_BRANCH}

pull_caars_env_docker:
	cd etc && ./build_caars_docker.sh only_pull_env
pull_caars_env_docker_travis:
	cd etc && ./build_caars_docker.sh only_pull_env ${TRAVIS_BRANCH}
pull_caars_docker:
	cd etc && ./build_caars_docker.sh only_pull

push_caars_env_docker:
	cd etc && ./build_caars_env.sh push_yes
push_caars_docker:
	cd etc && ./build_caars_docker.sh push_yes
push_caars_docker_travis:
	cd etc && ./build_caars_docker.sh push_yes ${TRAVIS_BRANCH}

##########################
# Some tests             #
##########################

# Basic test
test:
	cd example && bash Launch_CAARS.sh

##running test using a docker container by workflow in bistro
test_using_docker:
	cd example && bash Launch_CAARS.sh docker

clean_test:
	cd example && (rm -r working_dir/ output_dir/ || exit 0)

# test to check options
test_options:
	cd example && bash Launch_CAARS2.sh

##running test_options using a docker container by workflow in bistro
test_options_using_docker:
	cd example && bash Launch_CAARS2.sh docker

clean_test_options:
	cd example && (rm -r working2_dir/ output2_dir/ || exit 0)

# test to check options
test_wiki:
	mkdir test_wiki &&\
	cd test_wiki && \
	git clone https://github.com/CarineRey/caars.wiki.git . && \
	make tests

# test paper_pipline
test_paper_pipeline:
	dune build app/caars_paper_pipeline_app.exe
	mkdir -p test_paper_pipeline && \
	cd test_paper_pipeline && \
	cp ../_build/default/app/caars_paper_pipeline_app.exe caars_paper_pipeline_app && \
	./caars_paper_pipeline_app prepare-data --outdir dataset_preview --np 3 --preview && \
	./caars_paper_pipeline_app analysis --indir dataset_preview --outdir results --np 3


.PHONY: caars test_paper_pipeline test test_using_docker clean_test clean test_options test_options_using_docker clean_test_options build_caars_env_docker build_caars_docker push_caars_env_docker push_caars_docker

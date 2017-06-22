sudo apt-get update -qq
sudo apt-get install -y opam

opam init --comp="$OCAML_VERSION"
eval $(opam config env)

opam install -y bistro

make
make test

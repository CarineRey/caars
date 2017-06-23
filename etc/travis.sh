sudo apt-get update -qq
sudo apt-get install -y opam git

opam init --comp="$OCAML_VERSION"
eval $(opam config env)

opam install -y core=v0.9.1
opam install -y bistro

#git clone https://github.com/CarineRey/caars.git

make
make test

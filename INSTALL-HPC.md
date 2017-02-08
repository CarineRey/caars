# amalgam on HPC

**WORK IN PROGRESS**

## System-wide installation

Like standard locations, or distribution packages.

Some (rare) dependencies may be installed via distribution packages. Example, on Debian based distributions:

```sh
apt-get install python-setuptools python-qt4 python-scipy python-mysqldb python-lxml python-pip fasttree exonerate mafft

```

## Installation on non-standard locations

Here we describe a painfully, but complete, installation of amalgam and all it's dependencies (out of order). 

We will use [PSMN](http://www.ens-lyon.fr/PSMN/)'s case study as an example. PSMN use ```/applis/PSMN/``` as root filetree for all non-standard installations, on a NFS share on cluster's nodes. 

(TODO: user used in this document has read-write-execute rights on /applis/PSMN, explain differences between debian7/ & generic/)

### Python 2.7

With a site-packages on a non-standard location:

```sh
export PYTHONPATH="/applis/PSMN/debian7/python/2.7/site-packages"
export SITE="/applis/PSMN/debian7/python/2.7/site-packages" PREFIX="/applis/PSMN/debian7/python/2.7"

python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade pip
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade setuptools
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade PyQt4
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade scipy
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade MySQLdb
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade lxml
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade ete2
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade ete3
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade profileNJ
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade pandas
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade matplotlib
```

Cleanup:

```
export PREFIX=""
```

### Transdecoder >= 3.0.1

Get it from [https://github.com/TransDecoder/TransDecoder](https://github.com/TransDecoder/TransDecoder). It's written in Perl, but it build cd-hit, so there is a 'make' step.

```sh
mkdir -p /applis/PSMN/debian7/TransDecoder
cd /applis/PSMN/debian7/TransDecoder/
wget https://github.com/TransDecoder/TransDecoder/archive/v3.0.1.tar.gz
tar xvf v3.0.1.tar.gz
cd v3.0.1 && make
```

Add ```/applis/PSMN/debian7/TransDecoder/v3.0.1``` and ```/applis/PSMN/debian7/TransDecoder/v3.0.1/util/bin``` to ```PATH```.

### SRAToolKit >= 2.8.1-2

SRAToolKit is ditributed as binary. Get it from [https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software](https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=software)

```sh
mkdir -p /applis/PSMN/generic/SRAToolkit
cd /applis/PSMN/generic/SRAToolkitfasttree
wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.8.1-2/sratoolkit.2.8.1-2-centos_linux64.tar.gz
tar xzf sratoolkit.2.8.1-2-centos_linux64.tar.gz
```

Add ```/applis/PSMN/generic/SRAToolkit/sratoolkit.2.8.1-2-centos_linux64``` to ```PATH```.

### FastTree >= 2.1.7

FastTree is distributed as binary. Get it from [http://www.microbesonline.org/fasttree/](http://www.microbesonline.org/fasttree/). 

```sh
mkdir -p /applis/PSMN/generic/FastTree/2.1.9/
cd /applis/PSMN/generic/FastTree/2.1.9/
wget http://www.microbesonline.org/fasttree/FastTree
ln -s FastTree fasttree
```

Add ```/applis/PSMN/generic/FastTree/2.1.9/``` to ```PATH```.

### apytram >= 1.0

apytram is written in Python 2.7, but need following dependencies:

#### exonerate >= 2.2.0

exonerate is already packaged. See [System-wide installation](#system-wide-installation).

#### mafft >=7.1

mafft is already packaged. See [System-wide installation](#system-wide-installation).

#### blast+ >= 2.6

Blast+ is distributed as binary. Get it from [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi).

```
mkdir -p /applis/PSMN/generic/Blast/
wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz
tar zxf ncbi-blast-2.6.0+-x64-linux.tar.gz
mv blast-2.6.0+ 2.6.0+
```

Add ```/applis/PSMN/generic/Blast/2.6.0+/bin``` to ```PATH```.

#### Trinity >=2.3

Trinity can be downloaded on [github](https://github.com/trinityrnaseq/trinityrnaseq/releases). You'll need to install Java 1.8 (at least JRE 1.8u0) first.

##### Java >= 1.8

You can get Java from [java.com](http://javadl.oracle.com/webapps/download/).

```sh
mkdir -p /applis/PSMN/generic/java
cp jre-8u77-linux-x64.tar.gz /applis/PSMN/generic/java/
cd /applis/PSMN/generic/java/
tar --no-same-owner -xzvf jre-8u77-linux-x64.tar.gz
```

Add:
* ```/applis/PSMN/generic/java/jre1.8.0_77``` to ```JAVAHOME```,
* ```/applis/PSMN/generic/java/jre1.8.0_77/bin``` to ```PATH```.

TOBECONTINUED.

##### Bowtie >= 2

successfully tested with 2.2.9.

##### Building Trinity

#### Install apytram


### PhyloMerge (0.2 from 2017/01)
    * bpp >= 2.2.0 (Bio++)


### phyldog:
    * libPLL >= 1.0.2 sequential
    * boost 1.55 < . > 1.49 (waiting for PSMN's validation on 1.63)
    * bpp >= 2.2.0 (Bio++)

**WARNING**: phyldog, bpp, boost and mpi must be build with the same compiler.



### OCaml

Get latest OCaml from [http://ocaml.org/](http://ocaml.org/). Untar and build:

```sh
tar xvf ocaml-4.04.0.tar.gz
cd ocaml-4.04.0
./configure -prefix /applis/PSMN/debian7/OCaml/4.04.0
make world
make bootstrap
make opt
make install
```

Add:

* ```/applis/PSMN/debian7/OCaml/4.04.0/bin``` to ```PATH```, 
* ```/applis/PSMN/debian7/OCaml/4.04.0/lib``` to ```LD_LIBRARY_PATH```, 
* ```/applis/PSMN/debian7/OCaml/4.04.0/man``` to ```MANTPATH```. 

Use whatever-you-want environment tool :o) (as Module Environment, LMode, personnal rc files, etc.)

#### opam

* Install opam (OCaml Package Manager): 

See [https://opam.ocaml.org/](https://opam.ocaml.org/) for help.

```sh
wget https://raw.github.com/ocaml/opam/master/shell/opam_installer.sh
./opam_installer.sh /applis/PSMN/debian7/OCaml/4.04.0/bin
```

* Init opam

as standard user:

```sh
opam init
eval `opam config env`
```

* Install amalgam OCaml depencies:

as standard user:

```sh
opam pin add bistro --dev-repo
```

opam should automagically install bistro's dependencies:

* oasis
* solvuu-build
* ocamlgraph



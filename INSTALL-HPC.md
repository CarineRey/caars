# CAARS on HPC

**WORK IN PROGRESS**


## System-wide installation

Like standard locations, or distribution packages.

Some (rare) dependencies may be installed via distribution packages. Example, on Debian based distributions:

```sh
apt-get install python-setuptools python-qt4 python-scipy python-mysqldb python-lxml python-pip fasttree exonerate mafft openmpi-bin openmpi-checkpoint libopenmpi-dev

```

## Installation on non-standard locations

Here we describe a painfully, but complete, installation of CAARS and all it's dependencies (out of order). 

We will use [PSMN](http://www.ens-lyon.fr/PSMN/)'s case study as an example.

* PSMN use Debian 7 as it's base system,
* PSMN use ```/applis/PSMN/``` as root filetree for all non-standard installations, on a NFS share mounted on all cluster's nodes,
* PSMN discriminate between distributed binaries (/applis/PSMN/generic/) and compiled programs (/applis/PSMN/debian7/),
* PSMN use a dedicated user, with write permissions on ```/applis/PSMN/``` filetree, to install programs.

For each addition to environment variables (such as PATH, PYTHONPATH, etc.), use whatever-you-want environment tool :o) (such as Module Environment, LMode, personnal rc-files, etc.).

### Python 2.7

With a site-packages on PSMN's non-standard location:

```sh
export PYTHONPATH="/applis/PSMN/debian7/python/2.7/site-packages"
export LD_LIBRARY_PATH="/usr/lib/atlas-base:$LD_LIBRARY_PATH"
export SITE="/applis/PSMN/debian7/python/2.7/site-packages" PREFIX="/applis/PSMN/debian7/python/2.7"

python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade pip
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade setuptools
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade PyQt4
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade numpy
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade scipy
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade MySQLdb
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade lxml
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade ete2
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade ete3
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade profileNJ
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade pandas
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade matplotlib
python -m easy_install --prefix=$PREFIX --install-dir=$SITE --upgrade biopython
```

Cleanup:

```
export PREFIX=""
```

### Transdecoder >= 3.0.1

Get it from [github](https://github.com/TransDecoder/TransDecoder). It's written in Perl, but it build it's own 'cd-hit', so there is a 'make' step.

```sh
mkdir -p /applis/PSMN/debian7/TransDecoder
cd /applis/PSMN/debian7/TransDecoder/
wget https://github.com/TransDecoder/TransDecoder/archive/v3.0.1.tar.gz
tar xvf v3.0.1.tar.gz
cd v3.0.1 && make
```

Add ```/applis/PSMN/debian7/TransDecoder/v3.0.1``` and ```/applis/PSMN/debian7/TransDecoder/v3.0.1/util/bin``` to ```PATH```.

### FastTree >= 2.1.7

FastTree is distributed as binary. Get it from [microbesonline](http://www.microbesonline.org/fasttree/). 

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
mkdir -p /applis/PSMN/generic/Blast
cd /applis/PSMN/generic/Blast/
wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz
tar zxf ncbi-blast-2.6.0+-x64-linux.tar.gz
mv blast-2.6.0+ 2.6.0+
```

Add ```/applis/PSMN/generic/Blast/2.6.0+/bin``` to ```PATH```.

#### Trinity >=2.3

Trinity can be downloaded on [github](https://github.com/trinityrnaseq/trinityrnaseq/releases). 

You'll need to install both Java 1.8 (at least JRE 1.8u0) and Bowtie2 (successfully tested with 2.2.9) first.

##### Java >= 1.8

Get Java from [java.com](http://javadl.oracle.com/webapps/download/).

```sh
mkdir -p /applis/PSMN/generic/java
cp jre-8u77-linux-x64.tar.gz /applis/PSMN/generic/java/
cd /applis/PSMN/generic/java/
tar --no-same-owner -xzvf jre-8u77-linux-x64.tar.gz
```

Add:
* ```/applis/PSMN/generic/java/jre1.8.0_77``` to ```JAVAHOME```,
* ```/applis/PSMN/generic/java/jre1.8.0_77/bin``` to ```PATH```.

##### Bowtie >= 2

Bowtie2 is provided as binary and can be dowloaded on [sourceforge](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml).

```sh
mkdir -p /applis/PSMN/generic/Bowtie/2.2.9
cd /applis/PSMN/generic/Bowtie/2.2.9/
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.9/bowtie2-2.2.9-linux-x86_64.zip
unzip bowtie2-2.2.9-linux-x86_64.zip
mv bowtie2-2.2.9 x86_64
```

Add ```/applis/PSMN/generic/Bowtie/2.2.9/x86_64/``` to ```PATH```

##### Building Trinity

Trinity has to be builded "in-place", meaning there's no install, you have to build it where you want to install it.

```sh
mkdir -p /applis/PSMN/debian7/Trinity
cd /applis/PSMN/debian7/Trinity/
wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.3.2.tar.gz
tar zxf Trinity-v2.3.2.tar.gz
mv trinityrnaseq-Trinity-v2.3.2 2.3.2 && cd 2.3.2/
make
make plugins
```

Add ```/applis/PSMN/debian7/Trinity/2.3.2:/applis/PSMN/debian7/Trinity/2.3.2/trinity-plugins/```  to ```PATH```

#### Install apytram

The simpliest way is to get the latest development version of apytram, with git clone.

```sh
mkdir -p /applis/PSMN/generic/apytram
cd /applis/PSMN/generic/apytram/
git clone https://github.com/CarineRey/apytram.git dev
```

Add:
* ```/applis/PSMN/generic/apytram/dev``` to ```PATH```,
* ```/applis/PSMN/generic/apytram/dev/ApytramLib``` to ```PYTHONPATH```.

You can test if your apytram installation match all the requirements by using ```test_apytram_configuration.py```, or by running ```make test``` into apytram directory (/applis/PSMN/generic/apytram/dev/). See [apytram Wiki](https://github.com/CarineRey/apytram/wiki) for more informations.


### PhyloMerge (0.2 from 2017/01)

PhyloMerge has a dependency on Bio++ libraries.

#### Bio++

Get Bio++ installer from [univ-montp2.fr](http://biopp.univ-montp2.fr/wiki/index.php/Main_Page). Here's the version 2.2.0 installation example.

```sh
mkdir -p /applis/PSMN/debian7/Libs/bpp/2.2.0/gcc-4.7.2
mkdir -p ~/builds/biopp
cd ~/builds/biopp/
wget http://biopp.univ-montp2.fr/Download/bpp-setup.sh
```

Modify ```PATH_INSTALL``` in ```bpp-setup.sh``` file with your favorite text editor (MS PowerPoint©, for example).

```
PATH_INSTALL=/applis/PSMN/debian7/Libs/bpp/2.2.0/gcc-4.7.2
```

Then run ```./bpp-setup.sh```. It will download, build and install Bio++ libraries.

Add:
* ```/applis/PSMN/debian7/Libs/bpp/2.2.0/gcc-4.7.2/bin``` to ```PATH```,
* ```/applis/PSMN/debian7/Libs/bpp/2.2.0/gcc-4.7.2/include``` to ```INCLUDE``` and```CPATH```,
* ```/applis/PSMN/debian7/Libs/bpp/2.2.0/gcc-4.7.2/lib``` to ```LD_LIBRARY_PATH```, ```LD_RUN_PATH``` and```LIBRARY_PATH```.

#### Install PhyloMerge

Get PhyloMerge latest version from [github](https://github.com/boussau/phylomerge/).

```sh
mkdir -p /applis/PSMN/debian7/PhyloMerge/0.2/
cd ~/builds/
git clone https://github.com/boussau/phylomerge/ phylomerge
cd phylomerge/
g++ -s -pipe -o phylomerge PhyloMerge.cpp -I/applis/PSMN/debian7/Libs/bpp/2.2.0/gcc-4.7.2/include -L. -L/applis/PSMN/debian7/Libs/bpp/2.2.0/gcc-4.7.2/lib -O3 -fopenmp -std=c++0x /applis/PSMN/debian7/Libs/bpp/2.2.0/gcc-4.7.2/lib/libbpp-phyl.a /applis/PSMN/debian7/Libs/bpp/2.2.0/gcc-4.7.2/lib/libbpp-seq.a /applis/PSMN/debian7/Libs/bpp/2.2.0/gcc-4.7.2/lib/libbpp-core.a --static
```
After successfull build:

```sh
cp phylomerge /applis/PSMN/debian7/PhyloMerge/0.2/
```

Add ```/applis/PSMN/debian7/PhyloMerge/0.2/``` to ```PATH```.

### PHYLDOG 1.1.0

**WARNING**: PHYLDOG, bpp, boost and mpi must be build with the same compiler. On Debian 7, it's gcc 4.7.2.

#### OpenMPI >= 1.4.5

OpenMPI is already packaged. See [System-wide installation](#system-wide-installation).

#### Bio++ >= 2.2.0

Bio++ is already installed! See [Bio++](Bio++).

#### libPLL >= 1.0.2

The 'Phylogenetic Likelihood Library' can be downloaded from [libpll.org](http://www.libpll.org/). PHYLDOG works well with sequential version.

```sh
mkdir /applis/PSMN/generic/Libs/libpll/1.0.2/sse3
cd /applis/PSMN/generic/Libs/libpll/1.0.2/sse3/
wget http://www.libpll.org/Downloads/libpll-1.0.2-sse3-64.tar.gz
tar xzf libpll-1.0.2-sse3-64.tar.gz
mv libpll-1.0.2-sse3-64 seq
```

Add:
* ```/applis/PSMN/generic/Libs/libpll/1.0.2/sse3/seq/include``` to ```INCLUDE``` and```CPATH```,
* ```/applis/PSMN/generic/Libs/libpll/1.0.2/sse3/seq/lib``` to ```LD_LIBRARY_PATH```, ```LD_RUN_PATH``` and```LIBRARY_PATH```.

#### Boost > 1.49.0 & <= 1.55.0

Get Boost libraries form [boost.org](http://www.boost.org/). You will need an old version, latest know working is [1.55.0](https://sourceforge.net/projects/boost/files/boost/1.55.0/).


```sh
mkdir -p /applis/PSMN/debian7/Libs/Boost/1.55.0/gcc-4.7.2
cd ~/builds/
wget https://downloads.sourceforge.net/project/boost/boost/1.55.0/boost_1_55_0.tar.bz2
tar -zjvf boost_1_55_0.tar.bz2
cd boost_1_55_0/
./bootstrap.sh --prefix=/applis/PSMN/debian7/Libs/Boost/1.55.0/gcc-4.7.2 --with-libraries=all
./b2
./b2 install --prefix=/applis/PSMN/debian7/Libs/Boost/1.55.0/gcc-4.7.2
```

Add:
* ```/applis/PSMN/debian7/Libs/Boost/1.55.0/gcc-4.7.2``` to ```BOOST_ROOT```
* ```/applis/PSMN/debian7/Libs/Boost/1.55.0/gcc-4.7.2/include``` to ```INCLUDE``` and```CPATH```,
* ```/applis/PSMN/debian7/Libs/Boost/1.55.0/gcc-4.7.2lib``` to ```LD_LIBRARY_PATH```, ```LD_RUN_PATH``` and```LIBRARY_PATH```.


#### Install PHYLDOG

Get PHYLDOG from [github](https://github.com/Boussau/PHYLDOG).

```sh
mkdir -p /applis/PSMN/debian7/PHYLDOG/1.1.0/gcc-4.7.2/bin
cd ~/builds/
git clone https://github.com/Boussau/PHYLDOG phyldog
cd phyldog/
mkdir build && cd build/
cmake -DBUILD_STATIC=ON .. -DCMAKE_LIBRARY_PATH="/applis/PSMN/debian7/Libs/bpp/2.2.0/gcc-4.7.2/lib;/applis/PSMN/generic/Libs/libpll/1.0.2/sse3/seq/lib;/usr/lib/openmpi/lib" -DCMAKE_INCLUDE_PATH="/applis/PSMN/debian7/Libs/bpp/2.2.0/gcc-4.7.2/include;/applis/PSMN/generic/Libs/libpll/1.0.2/sse3/seq/include;/usr/lib/openmpi/include" -DBOOST_LIBRARYDIR="/applis/PSMN/debian7/Libs/Boost/1.55.0/gcc-4.7.2/lib" -DBOOST_ROOT="/applis/PSMN/debian7/Libs/Boost/1.55.0/gcc-4.7.2/" -DCMAKE_INSTALL_PREFIX=/applis/PSMN/debian7/PHYLDOG/1.1.0/gcc-4.7.2
make
make install
```

Add ```/applis/PSMN/debian7/PHYLDOG/1.1.0/gcc-4.7.2/bin``` to ```PATH```.

### OCaml >= 4.03.0

Get latest OCaml from [ocaml.org](http://ocaml.org/). Untar and build:

```sh
tar xvf ocaml-4.04.0.tar.gz
cd ocaml-4.04.0/
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


#### opam

opam is OCaml Package Manager. See [ocaml.org](https://opam.ocaml.org/) for help.

```sh
wget https://raw.github.com/ocaml/opam/master/shell/opam_installer.sh
./opam_installer.sh /applis/PSMN/debian7/OCaml/4.04.0/bin
```

* Init opam as standard user:

```sh
opam init
eval `opam config env`
```

* Install CAARS's OCaml depencies as standard user:

```sh
opam pin add bistro --dev-repo
```

opam should automagically install bistro's dependencies:

* oasis
* solvuu-build
* ocamlgraph

### CAARS developement version

Finally (at last), install CAARS from [github](https://github.com/CarineRey/caars).

```sh
mkdir -p /applis/PSMN/debian7/caars
cd /applis/PSMN/debian7/caars/
git clone https://github.com/carinerey/caars dev
cd dev/
make
```

Add:
* ```/applis/PSMN/debian7/caars/utils/bin:/applis/PSMN/debian7/caars``` to ```PATH```,
* ```/applis/PSMN/debian7/caars/lib``` to ```PYTHONPATH```.

Et voilà! You're done.

TO DO:

* seqtk
* biopython
* cd-hit
* remove sratoolkit

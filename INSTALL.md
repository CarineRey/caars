# CAARS 

## Dependencies tree

* Transdecoder >= 3.0.1

* FastTree >= 2.1.7 (Warning: The executable must be fasttree and not FastTree)

* Python 2.7 (with pip and setuptools)
    * PyQt4
    * SciPy
    * MySQLdb
    * lxml
    * ete2
    * ete3
    * profileNJ
    * pandas
    * Biopython

* Trinity >=2.3
    * Java >= 1.8 (OpenJRE works)
    * Bowtie >= 2 (tested with 2.2.9)

* phyldog >= 1.1.0
    * libPLL >= 1.0.2 sequential
    * boost from 1.49 to 1.55 (versions >1.55 won't do for now)
    * bpp >= 2.2.0 (Bio++)

* PhyloMerge (0.2 from 2017/01)
    * bpp >= 2.2.0 (Bio++)

* apytram >= 1.0
    * exonerate >= 2.2.0
    * mafft >=7.1
    * blast+ >= 2.6
    * Seqtk >= 1.2
    * python = 2.7
    * cd-hit >= 4.6.6
    * Trinity >=2.3

* OCaml >= 4.03.0
    * bistro = 0.3.0


## dependencies installation via distribution packages

Like standard locations, or distribution packages.

Some (rare) dependencies may be installed via distribution packages. Example, on Debian based distributions:

```sh
apt-get install python-setuptools python-qt4 python-scipy python-mysqldb python-lxml python-pip fasttree exonerate mafft

```
You may need root rights (use for instance ```sudo```)

## Installation on non-standard locations

Here we describe a painfull, but complete, installation of CAARS and all its dependencies (out of order).

### Python 2.7

With a site-packages on a non-standard location:

```sh
python -m easy_install --upgrade pip
python -m easy_install --upgrade setuptools
python -m easy_install --upgrade PyQt4
python -m easy_install --upgrade scipy
python -m easy_install --upgrade MySQLdb
python -m easy_install --upgrade lxml
python -m easy_install --upgrade ete2
python -m easy_install --upgrade ete3
python -m easy_install --upgrade profileNJ
python -m easy_install --upgrade pandas
python -m easy_install --upgrade matplotlib
python -m easy_install --upgrade biopython
```

### Transdecoder >= 3.0.1

Get it from [https://github.com/TransDecoder/TransDecoder](https://github.com/TransDecoder/TransDecoder). It's written in Perl, but it builds cd-hit.

``` sh
mkdir -p /home/user/bin/TransDecoder
cd /home/user/bin/TransDecoder/
wget https://github.com/TransDecoder/TransDecoder/archive/v3.0.1.tar.gz
tar xvf v3.0.1.tar.gz
cd v3.0.1 && make
```

Add ```/home/user/bin/TransDecoder/v3.0.1``` and ```/home/user/bin/TransDecoder/v3.0.1/util/bin``` to ```PATH```.

### FastTree >= 2.1.7

FastTree is distributed as binary. Get it from [http://www.microbesonline.org/fasttree/](http://www.microbesonline.org/fasttree/). 

```sh
mkdir -p /home/user/bin/FastTree/2.1.9/
cd /home/user/bin/FastTree/2.1.9/
wget http://www.microbesonline.org/fasttree/FastTree
ln -s FastTree fasttree
```

Add ```/home/user/bin/FastTree/2.1.9/``` to ```PATH```.

### apytram >= 1.0

apytram is written in Python 2.7, but need following dependencies:


#### exonerate >= 2.2.0

exonerate is already packaged. See [Dependencies installation via distribution packages](#dependencies-installation-via-distribution-packages).

#### mafft >=7.1

mafft is already packaged. See [Dependencies installation via distribution packages](#dependencies-installation-via-distribution-packages).

#### blast+ >= 2.6

Blast+ is distributed as binary, which can be downloaded at the [NCBI](https://blast.ncbi.nlm.nih.gov/Blast.cgi).

```
mkdir -p /home/user/bin/Blast/
wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/2.6.0/ncbi-blast-2.6.0+-x64-linux.tar.gz
tar zxf ncbi-blast-2.6.0+-x64-linux.tar.gz
mv blast-2.6.0+ 2.6.0+
```

Add ```/home/user/bin/Blast/2.6.0+/bin``` to ```PATH```.

#### Seqtk  >= 1.2

Toolkit for processing sequences in FASTA/Q formats [https://github.com/lh3/seqtk](https://github.com/lh3/seqtk)

``` sh
mkdir -p /home/user/bin/seqtk
git clone https://github.com/lh3/seqtk.git  /home/user/bin/seqtk
cd  /home/user/bin/seqtk
make
```

Add ```/home/user/bin/seqtk``` to ```PATH```.

``` sh
export PATH=/home/user/bin/seqtk:$PATH
```


#### CD-HIT >= 4.6.6


CD-HIT is a very widely used program for clustering and comparing protein or nucleotide sequences ([http://weizhongli-lab.org/cd-hit/](http://weizhongli-lab.org/cd-hit/)).


``` sh
 git clone https://github.com/weizhongli/cdhit.git /home/user/bin/cdhit/
cd /home/user/bin/cdhit/cd-hit-auxtools
make
cd /home/user/bin/cdhit
make
make install
```

Add ```/home/user/bin/cd-hit-auxtools``` to ```PATH```.

``` sh
export PATH=/home/user/bin/cd-hit-auxtools:$PATH
```



#### Install apytram

```sh
git clone https://github.com/carinerey/apytram /home/user/bin/apytram
```

Add:
*  ```/home/user/bin/apytram``` to ```PATH```.
*  ```/home/user/bin/apytram/ApytramLib``` to ```PYTHONPATH```.

You can test if your apytram installation match all the requirements by using ```test_apytram_configuration.py```,
or by running ```make test``` into apytram directory (/home/user/bin/apytram).
See [apytram Wiki](https://github.com/CarineRey/apytram/wiki) for more informations.




#### Trinity >=2.3

Trinity can be downloaded on [github](https://github.com/trinityrnaseq/trinityrnaseq/releases). You'll need to install first:

* Java 1.8 (at least JRE 1.8u0)
* Bowtie >= 2

##### Java >= 1.8

You can get Java from [java.com](http://javadl.oracle.com/webapps/download/).

```sh
mkdir -p /home/user/bin/java
cp jre-8u77-linux-x64.tar.gz /home/user/bin/java/
cd /home/user/bin/java/
tar --no-same-owner -xzvf jre-8u77-linux-x64.tar.gz
```

Add:
* ```/home/user/bin/java/jre1.8.0_77``` to ```JAVAHOME```,
* ```/home/user/bin/java/jre1.8.0_77/bin``` to ```PATH```.

##### Bowtie >= 2
successfully tested with 2.2.9.

```sh
cd /home/user/bin/
wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.0/bowtie2-2.3.0-source.zip
unzip bowtie2-2.3.0-source.zip
cd bowtie2-2.3.0/
make
make install
```


### PhyloMerge (0.2 from 2017/01)
    * bpp >= 2.2.0 (Bio++)


#### Bpp >= 2.2.0

Bio++ Program Suite  ([https://github.com/BioPP/bppsuite](https://github.com/BioPP/bppsuite))

``` sh
echo "deb http://biopp.univ-montp2.fr/repos/apt/ Trusty main" >> /etc/apt/sources.list;\
    wget http://biopp.univ-montp2.fr/repos/apt/conf/biopp.gpg.key
apt-key add biopp.gpg.key
apt-get update
apt-get install -qy libbpp-phyl-dev

```


#### Phylomerge 


``` sh
git clone https://github.com/boussau/phylomerge/ /home/user/bin/phylomerge
cd  /home/user/bin/phylomerge
bash ./compile
```


Add ```/home/user/bin/phylomerge``` to ```PATH```.

``` sh
export PATH=/home/user/bin/phylomerge:$PATH
```



### phyldog:
    * libPLL >= 1.0.2 sequential
    * boost 1.55 < . > 1.49 (waiting for PSMN's validation on 1.63)
    * bpp >= 2.2.0 (Bio++) (see above)

**WARNING**: phyldog, bpp, boost and mpi must be build with the same compiler.


#### libPLL >= 1.0.2

``` sh
mkdir -p /home/user/lib/
cd /home/user/lib/
wget http://www.libpll.org/Downloads/libpll-1.0.2-sse3-64.tar.gz
tar xzf libpll-1.0.2-sse3-64.tar.gz && rm libpll-1.0.2-sse3-64.tar.gz
cd /home/user/lib/libpll-1.0.2-sse3-64/
cp libpll-sse3* /usr/local/lib
cp -r include/pll /usr/local/include/
```

#### Phyldog

``` sh
cd /home/user/bin/
git clone https://github.com/Boussau/PHYLDOG phyldog

mkdir -p /home/user/bin/phyldog/build
cd /home/user/bin/phyldog/build


export LD_LIBRARY_PATH=/usr/local/lib/:$LD_LIBRARY_PATH  #(where you store libpll)

cmake ..
make
make install
```

### OCaml

Get latest OCaml from [http://ocaml.org/](http://ocaml.org/). Untar and build:

```sh
tar xvf ocaml-4.04.0.tar.gz
cd ocaml-4.04.0
./configure -prefix /home/user/binOCaml/4.04.0
make world
make bootstrap
make opt
su make install
```

Add:

* ```/home/user/binOCaml/4.04.0/bin``` to ```PATH```, 
* ```/home/user/binOCaml/4.04.0/lib``` to ```LD_LIBRARY_PATH```, 
* ```/home/user/binOCaml/4.04.0/man``` to ```MANTPATH```. 

Use whatever-you-want environment tool :o) (as Module Environment, LMode, personnal rc files, etc.)

#### opam

* Install opam (OCaml Package Manager): 

See [https://opam.ocaml.org/](https://opam.ocaml.org/) for help.

```sh
wget https://raw.github.com/ocaml/opam/master/shell/opam_installer.sh
su ./opam_installer.sh /home/user/bin/OCaml/4.04.0/bin
```

* Init opam

as user:

```sh
opam init
eval `opam config env`
```

* Install CAARS OCaml depencies:

as user:

```sh
opam install -y bistro=0.3.0
```

### CAARS

Finally (at last), install CAARS from [github](https://github.com/CarineRey/caars).


```sh
git clone https://github.com/carinerey/caars /home/user/bin/caars/
cd  /home/user/bin/caars/
make
```

Add ```/home/user/bin/caars/``` to ```PATH```.
Add ```/home/user/bin/caars/utils/bin/``` to ```PATH```.

``` sh
export PATH=/home/user/bin/caars/:$PATH
```

To test that all is okay, you can run `make test`
``` sh
cd  /home/user/bin/caars/
make test
```

FROM debian:buster-slim
MAINTAINER Carine Rey carine.rey@ens-lyon.org

RUN apt-get clean && apt-get update

RUN apt-get install --no-install-recommends  -qy git \
                        wget \
                        cmake \
                        build-essential \
                        mafft \
                        exonerate \
                        python-pip \
                        python-dev \
                        python-setuptools \
                        libidn11 \
# to install trinity
                        zlib1g-dev \
                        ncurses-dev \
                        rsync \
                        samtools \
                        jellyfish \
                        salmon \
#To install Transdecoder
                        libany-uri-escape-perl \
# to install bowtie2
                        libtbb-dev \
                        unzip \
# to user
                        vim \
                        htop \
                        xvfb \
                        gosu
RUN apt-get install --no-install-recommends  -qy \
# to mafft
                        ruby \
#to install ete3
                        python-setuptools \
                        python-qt4 \
                        python-tk \
                        python-scipy \
                        python-mysqldb \
                        python-lxml \
#To install opam
                        bzip2 \
                        curl \
                        make \
                        m4 \
                        gcc \
                        pkg-config \
                        libncurses5-dev \
                        libgdbm-dev \
# to install generax
                        flex \
                        bison \
                        libgmp3-dev \
                        g++ \
                        libboost-all-dev \
                        ssh \ 
#mpi needs something to communicate (rsh, ssh other)
# to install fasttree
                        fasttree \
                        sudo \
                        aspcud \
                        less

RUN pip install wheel numpy pandas matplotlib biopython

#### install
RUN pip install ete3==3.0.0b35

### install openjdk-11-jre
RUN echo "deb http://httpredir.debian.org/debian stretch-backports main" >> /etc/apt/sources.list && mkdir /usr/share/man/man1 && \
    apt-get update && apt-get install --no-install-recommends -t stretch-backports -qy ca-certificates-java openjdk-11-jre-headless 


###  install bowtie2 2.3.5.1
ENV bowtie_version 2.3.5.1
WORKDIR /opt
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/"$bowtie_version"/bowtie2-"$bowtie_version"-source.zip &&\
    unzip bowtie2-"$bowtie_version"-source.zip && cd bowtie2-"$bowtie_version"/ && make && make install


###  install ncbi-blast+=2.9.0
ENV blast_version 2.9.0
WORKDIR /opt
RUN wget ftp://ftp.ncbi.nih.gov/blast/executables/blast+/"$blast_version"/ncbi-blast-"$blast_version"+-x64-linux.tar.gz &&\
    tar zxf ncbi-blast-"$blast_version"+-x64-linux.tar.gz
ENV PATH /opt/ncbi-blast-"$blast_version"+/bin/:$PATH


### install Trinity
ENV trinity_version v2.9.1
WORKDIR /opt

RUN wget https://github.com/trinityrnaseq/trinityrnaseq/releases/download/"$trinity_version"/trinityrnaseq-"$trinity_version".FULL.tar.gz && \
     tar zxf trinityrnaseq-"$trinity_version".FULL.tar.gz && rm trinityrnaseq-"$trinity_version".FULL.tar.gz &&\
     cd /opt/trinityrnaseq-"$trinity_version"/ && make && make install
ENV PATH $PATH:/opt/trinityrnaseq-"$trinity_version"/:/opt/trinityrnaseq-"$trinity_version"/trinity-plugins/


### install seqtk
WORKDIR /opt/seqtk
RUN git clone https://github.com/lh3/seqtk.git  /opt/seqtk/ && make
ENV PATH /opt/seqtk:$PATH


### install cd-hit
RUN git clone https://github.com/weizhongli/cdhit.git /opt/cdhit/
RUN cd /opt/cdhit/cd-hit-auxtools &&  make && \
    cd /opt/cdhit &&  make  && make install
ENV PATH /opt/cdhit/cd-hit-auxtools:$PATH


#### install TransDecoder
ENV transdecoder_version TransDecoder-v5.5.0
WORKDIR /opt

RUN wget https://github.com/TransDecoder/TransDecoder/archive/"$transdecoder_version".zip &&\
    unzip "$transdecoder_version".zip  && rm "$transdecoder_version".zip && \
    ls /opt/ && cd /opt/TransDecoder-"$transdecoder_version"/ && make && ls /opt
ENV PATH /opt/TransDecoder-"$transdecoder_version"/:$PATH

## install R and ggplot2 for TransDecoder
#RUN echo 'install.packages(c("ggplot2"), repos="http://cran.us.r-project.org", dependencies=TRUE); \
#    source("https://bioconductor.org/biocLite.R");\
#    biocLite("seqLogo")' > /tmp/install.R && Rscript /tmp/install.R && rm /tmp/install.R

### install apytram
WORKDIR /opt/apytram
RUN git clone https://github.com/carinerey/apytram /opt/apytram &&\
    git pull && git checkout a3cb5da
ENV PATH /opt/apytram:$PATH

### Install muscle
WORKDIR /usr/local/bin
RUN wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz && \
    tar xvzf muscle3.8.31_i86linux64.tar.gz && rm muscle3.8.31_i86linux64.tar.gz && mv muscle3.8.31_i86linux64 muscle


##### install generax
WORKDIR  /opt/
RUN git clone --recursive https://github.com/BenoitMorel/GeneRax && cd GeneRax && ./install.sh
ENV PATH /opt/GeneRax/build/bin/:$PATH

#### install bpp

ENV bpp_dir=/usr/local/lib

WORKDIR $bpp_dir/sources_bpp
RUN git clone  https://github.com/BioPP/bpp-core &&\
    cd  $bpp_dir/sources_bpp/bpp-core &&\
    git checkout 1020f99f55fc6e122ff5d70b5ef4bae272307951 &&\
    cmake . &&\
    (make -j 4 || make) &&\
    make install
    
WORKDIR $bpp_dir/sources_bpp
RUN git clone https://github.com/BioPP/bpp-seq &&\
    cd $bpp_dir/sources_bpp/bpp-seq &&\
    git checkout 42f56bbd4802c08c30d92c521c58eef49793d4bd &&\
    cmake  . &&\
    (make -j 4 || make) &&\
    make install

WORKDIR $bpp_dir/sources_bpp
RUN git clone https://github.com/BioPP/bpp-phyl &&\
    cd  $bpp_dir/sources_bpp/bpp-phyl &&\
    git checkout ee223fb9b937c6dc7aab457547c2912c1bf434c0 &&\
    cmake - . &&\
    (make -j 4 || make) &&\
    make install &&\
    rm -r $bpp_dir/sources_bpp


##### install phylomerge
RUN git clone https://github.com/boussau/phylomerge/ /opt/phylomerge &&\
    cd /opt/phylomerge && bash ./compile
ENV PATH /opt/phylomerge:$PATH


### install opam
WORKDIR /home/caars_user
ENV HOME /home/caars_user
RUN apt install -y opam

RUN opam init --comp 4.09.0 --disable-sandboxing
ENV PATH /usr/local/bin/.opam/4.09.0/bin:$PATH
RUN echo ". $HOME/.opam/opam-init/init.sh > /dev/null 2> /dev/null || true" >> ~/.profile

ENV PATH $HOME/.opam/4.09.0/bin/:$PATH
ENV PATH $HOME/.opam/4.09.0/bin/:$PATH
ENV CAML_LD_LIBRARY_PATH $HOME/.opam/4.09.0/lib/stublibs
ENV OCAML_TOPLEVEL_PATH $HOME/.opam/4.09.0/lib/toplevel
ENV MANPATH $HOME/.opam/4.09.0/man:$MANPATH
ENV PERL5LIB $HOME/.opam/4.09.0/lib/perl5:$PERL5LIB

##### install bistro
RUN opam pin add -y bistro --dev-repo
RUN opam install ocamlify

##### install local dependencies
WORKDIR /opt/caars_py_lib
ENV PYTHONPATH  /opt/caars_py_lib:$PYTHONPATH
COPY utils/lib/*.py /opt/caars_py_lib/
COPY utils/bin/* /usr/bin/

WORKDIR /opt/

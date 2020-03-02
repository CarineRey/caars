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
# to install trinity
                        zlib1g-dev \
                        ncurses-dev \
                        rsync \
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


RUN pip install numpy
RUN pip install pandas
RUN pip install matplotlib
RUN pip install biopython

### install openjdk-11-jre
RUN echo "deb http://httpredir.debian.org/debian stretch-backports main" >> /etc/apt/sources.list && mkdir /usr/share/man/man1 && \
    apt-get update && apt-get install --no-install-recommends -t stretch-backports -qy ca-certificates-java openjdk-11-jre-headless 


###  install bowtie2 2.3.3.1
ENV bowtie_version 2.3.3.1
WORKDIR /opt
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/"$bowtie_version"/bowtie2-"$bowtie_version"-source.zip &&\
    unzip bowtie2-"$bowtie_version"-source.zip && cd bowtie2-"$bowtie_version"/ && make && make install


###  install ncbi-blast+=2.6.0
ENV blast_version 2.6.0
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
ENV transdecoder_version 3.0.1
WORKDIR /opt
RUN wget https://github.com/TransDecoder/TransDecoder/archive/v"$transdecoder_version".zip &&\
    unzip v"$transdecoder_version".zip  && rm v"$transdecoder_version".zip && \
    ls /opt/ && cd /opt/TransDecoder-"$transdecoder_version"/ && make && ls /opt
ENV PATH /opt/TransDecoder-"$transdecoder_version"/:$PATH

## install R and ggplot2 for TransDecoder
#RUN echo 'install.packages(c("ggplot2"), repos="http://cran.us.r-project.org", dependencies=TRUE); \
#    source("https://bioconductor.org/biocLite.R");\
#    biocLite("seqLogo")' > /tmp/install.R && Rscript /tmp/install.R && rm /tmp/install.R

### install apytram
RUN git clone https://github.com/carinerey/apytram /opt/apytram
ENV PATH /opt/apytram:$PATH

#### install ete2
RUN apt-get install --no-install-recommends -qy unzip \
                        ruby \
#to install ete2
                        python-setuptools \
                        python-qt4 \
                        python-tk \
                        python-scipy \
                        python-mysqldb \
                        python-lxml \
#To install Transdecoder
                        libany-uri-escape-perl \
#To install opam
                        bzip2 \
                        curl \
                        make \
                        m4 \
                        gcc \
                        pkg-config \
                        libncurses5-dev \
                        libgdbm-dev \
# to install phyldog
                        cmake \
                        g++ \
                        libboost-all-dev \
                        ssh \ 
#mpi needs something to communicate (rsh, ssh other)
# to install fasttree
                        fasttree \
                        sudo \
                        aspcud \
                        less

RUN pip install ete2

### Install muscle
WORKDIR /usr/local/bin
RUN wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz && \
    tar xvzf muscle3.8.31_i86linux64.tar.gz && rm muscle3.8.31_i86linux64.tar.gz && mv muscle3.8.31_i86linux64 muscle

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

### install opam
WORKDIR /home/caars_user
ENV HOME /home/caars_user
RUN apt install -y opam
ENV PATH /usr/local/bin/.opam/4.09.0/bin:$PATH

RUN opam init --comp 4.09.0 --disable-sandboxing
RUN echo ". $HOME/.opam/opam-init/init.sh > /dev/null 2> /dev/null || true" >> ~/.profile

ENV PATH $HOME/.opam/4.09.0/bin/:$PATH
ENV CAML_LD_LIBRARY_PATH $HOME/.opam/4.09.0/lib/stublibs
ENV OCAML_TOPLEVEL_PATH $HOME/.opam/4.09.0/lib/toplevel
ENV MANPATH $HOME/.opam/4.09.0/man:$MANPATH
ENV PERL5LIB $HOME/.opam/4.09.0/lib/perl5:$PERL5LIB

##### install bistro
RUN opam pin add -y bistro --dev-repo

#### install profileNJ
RUN pip install --upgrade pip
RUN pip install ete3==3.0.0b35
RUN pip install profileNJ

##### install phylomerge
RUN git clone https://github.com/boussau/phylomerge/ /opt/phylomerge &&\
    cd /opt/phylomerge && bash ./compile
ENV PATH /opt/phylomerge:$PATH

## update apytram
WORKDIR /opt/apytram
RUN git pull && git checkout a3cb5da

##### install phyldog
WORKDIR  /opt/phyldog/

# for libpll
RUN apt-get clean && apt-get update && apt-get install --no-install-recommends -qy flex bison autotools-dev autoconf libtool automake

######### pll
COPY libpll-1.0.2-sse3-64.tar.gz /opt/libpll-1.0.2-sse3-64.tar.gz
WORKDIR /opt/
RUN tar xzf libpll-1.0.2-sse3-64.tar.gz && rm libpll-1.0.2-sse3-64.tar.gz
RUN cd /opt/libpll-1.0.2-sse3-64/ && cp libpll-sse3* /usr/local/lib && cp -r include/pll /usr/local/include/

######### libpll2
WORKDIR  /opt/phyldog/
RUN git clone --branch faster_reconciliation  https://github.com/BenoitMorel/PHYLDOG.git /opt/phyldog/ &&\
    git checkout bea224f1d46ed17474627f33774f2a9798ef38d9
RUN sed -e "s|git@github.com:ddarriba/pll-modules.git| --branch dev https://github.com/ddarriba/pll-modules.git|" -i ./build_libpll2.sh && ./build_libpll2.sh
RUN mv /opt/phyldog/libpll2_dependencies/include/pll /opt/phyldog/libpll2_dependencies/include/pllmodules
RUN mkdir -p build && cd /opt/phyldog/build && cmake  -DCMAKE_LIBRARY_PATH=/opt/phyldog/libpll2_dependencies/libs -DCMAKE_INCLUDE_PATH=/opt/phyldog/libpll2_dependencies/include ..   && (make -j 4 install || make install)
ENV PATH /opt/phyldog/build/bin/:$PATH

##### install local dependencies
WORKDIR /opt/caars_py_lib
ENV PYTHONPATH  /opt/caars_py_lib:$PYTHONPATH
COPY utils/lib/*.py /opt/caars_py_lib/
COPY utils/bin/* /usr/bin/

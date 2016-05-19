FROM carinerey/apytram
MAINTAINER Carine Rey carine.rey@ens-lyon.fr

#### install ete2
RUN apt-get install -qy python-setuptools \
                        python-qt4 \
                        python-scipy \
                        python-mysqldb \
                        python-lxml \
                        python-pip

RUN pip install ete2

#### install bpp
WORKDIR /opt
RUN echo "deb http://biopp.univ-montp2.fr/repos/apt/ Trusty main" >> /etc/apt/sources.list;\
    wget http://biopp.univ-montp2.fr/repos/apt/conf/biopp.gpg.key; \
    apt-key add biopp.gpg.key; \
	apt-get update ; \
    apt-get install -qy libbpp-phyl-dev 

##### install phyldog
RUN apt-get install -y cmake \
                       g++-4.9

RUN git clone git://pbil-dev.univ-lyon1.fr/phyldog /opt/phyldog && \
    mkdir /opt/phyldog/build && \
    apt-get install -y libboost-all-dev

WORKDIR /root/usr/bin
RUN wget http://www.libpll.org/Downloads/libpll-1.0.11.tar.gz && \
    tar xvzf libpll-1.0.11.tar.gz

WORKDIR /root/usr/bin/libpll-1.0.11
ENV CXX /usr/bin/g++-4.9
RUN ./configure && make install 

WORKDIR /opt/phyldog/build
RUN cmake .. -DCMAKE_CXX_COMPILER=/usr/bin/g++-4.9 && make install 

##### install fastree
RUN apt-get install -qy fasttree  # pranck

##### install opam 
RUN apt-get update && apt-get install -y \
    aspcud bzip2 curl wget git  curl make m4 unzip gcc \
    pkg-config libncurses5-dev libgdbm-dev \
    sudo
 
RUN wget https://raw.github.com/ocaml/opam/master/shell/opam_installer.sh -O - | sh -s /root/usr/bin 4.02.3
ENV PATH /root/usr/bin/.opam/4.02.3/bin:$PATH
RUN /root/usr/bin/opam init --comp 4.02.3

RUN echo ". /root/usr/bin/.opam/opam-init/init.sh > /dev/null 2> /dev/null || true" >> ~/.profile
RUN echo "export PATH=/root/usr/bin/:$PATH" >> ~/.profile


ENV PATH /root/usr/bin/:/root/.opam/4.02.3/bin/:$PATH
ENV CAML_LD_LIBRARY_PATH /root/.opam/4.02.3/lib/stublibs
ENV OCAML_TOPLEVEL_PATH /root/.opam/4.02.3/lib/toplevel
ENV MANPATH /root/.opam/4.02.3/man
ENV PERL5LIB /root/.opam/4.02.3/lib/perl5

#### install bistro
RUN opam pin add sexplib 113.00.00
RUN opam pin add bistro --dev-repo

##### install phylomerge
RUN git clone https://github.com/boussau/phylomerge/ /opt/phylomerge
WORKDIR /opt/phylomerge
RUN bash ./compile
ENV PATH /opt/phylomerge:$PATH


#### install amalgam

RUN git clone https://github.com/carinerey/amalgam /root/mydisk/amalgam/

WORKDIR /root/mydisk/amalgam/Bistro/
RUN  bash compile_commands.sh    
ENV PATH /root/mydisk/amalgam/Bistro:$PATH

#### install sratoolkit
WORKDIR /opt/
RUN wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-ubuntu64.tar.gz"
RUN tar -zxf sratoolkit.2.5.2-ubuntu64.tar.gz
RUN rm sratoolkit.2.5.2-ubuntu64.tar.gz
ENV  PATH /opt/sratoolkit.2.5.2-ubuntu64/bin/:$PATH

#### build a workingd irectory
RUN mkdir -p /root/mydisk/amalgam/data
WORKDIR /root/mydisk/amalgam/data


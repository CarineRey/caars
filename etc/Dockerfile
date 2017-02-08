FROM carinerey/apytram
MAINTAINER Carine Rey carine.rey@ens-lyon.fr

RUN apt-get clean && apt-get update

#### install ete2
RUN apt-get install --no-install-recommends -qy unzip \
#to install ete2
                        python-setuptools \
                        python-qt4 \
                        python-scipy \
                        python-mysqldb \
                        python-lxml \
                        python-pip \
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
                        g++-4.9 \
                        libboost-all-dev


RUN pip install ete2

##### install fasttree
RUN apt-get install -qy fasttree

#### install sratoolkit
WORKDIR /opt/
RUN wget "http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.5.2/sratoolkit.2.5.2-ubuntu64.tar.gz"
RUN tar zxf sratoolkit.2.5.2-ubuntu64.tar.gz
RUN rm sratoolkit.2.5.2-ubuntu64.tar.gz
ENV  PATH /opt/sratoolkit.2.5.2-ubuntu64/bin/:$PATH

#### install TransDecoder
WORKDIR /opt
RUN wget https://github.com/TransDecoder/TransDecoder/archive/v3.0.0.zip
RUN unzip v3.0.0.zip  && rm v3.0.0.zip
WORKDIR /opt/TransDecoder-3.0.0/
RUN make
ENV PATH /opt/TransDecoder-3.0.0/:$PATH

#### install opam
WORKDIR /root/
RUN wget https://raw.github.com/ocaml/opam/master/shell/opam_installer.sh -O - | sh -s /root/usr/bin 4.03.0
ENV PATH /root/usr/bin/.opam/4.03.0/bin:$PATH
RUN /root/usr/bin/opam init --comp 4.03.0

RUN echo ". /root/usr/bin/.opam/opam-init/init.sh > /dev/null 2> /dev/null || true" >> ~/.profile
RUN echo "export PATH=/root/usr/bin/:$PATH" >> ~/.profile

ENV PATH /root/usr/bin/:/root/.opam/4.03.0/bin/:$PATH
ENV CAML_LD_LIBRARY_PATH /root/.opam/4.03.0/lib/stublibs
ENV OCAML_TOPLEVEL_PATH /root/.opam/4.03.0/lib/toplevel
ENV MANPATH /root/.opam/4.03.0/man
ENV PERL5LIB /root/.opam/4.03.0/lib/perl5

##### install bistro
#RUN opam pin add -y oasis git://github.com/ocaml/oasis.git#opam/testing
#RUN opam pin add -y solvuu-build https://github.com/solvuu/solvuu-build.git#45dfffef8587d45ec9166615418749372ff4a0e8
#RUN opam install ocamlgraph
RUN opam pin add bistro --dev-repo


#### install bpp
WORKDIR /opt
RUN echo "deb http://biopp.univ-montp2.fr/repos/apt/ Trusty main" >> /etc/apt/sources.list;\
    wget http://biopp.univ-montp2.fr/repos/apt/conf/biopp.gpg.key &&\
    apt-key add biopp.gpg.key && \
	apt-get update && \
    apt-get install -qy libbpp-phyl-dev


##### install phylomerge
RUN git clone https://github.com/boussau/phylomerge/ /opt/phylomerge
WORKDIR /opt/phylomerge
RUN bash ./compile
ENV PATH /opt/phylomerge:$PATH

#### install profileNJ
RUN pip install --upgrade pip
RUN pip install ete3
RUN pip install profileNJ


#### install amalgam
RUN git clone https://github.com/carinerey/amalgam /usr/local/amalgam/
WORKDIR /usr/local/amalgam/
RUN make
ENV PATH /usr/local/amalgam/utils/bin:/usr/local/amalgam:$PATH
ENV PYTHONPATH /usr/local/amalgam/utils/lib



##### install phyldog

######### libpll-1.0.11
##  WORKDIR /root/usr/bin
##  RUN wget http://www.libpll.org/Downloads/libpll-1.0.11.tar.gz && \
##      tar xzf libpll-1.0.11.tar.gz
##  
##  WORKDIR /root/usr/bin/libpll-1.0.11
##  ENV CXX /usr/bin/g++-4.9
##  RUN ./configure && make install
##  
##  WORKDIR /opt/phyldog/build
##  RUN cmake .. -DCMAKE_CXX_COMPILER=/usr/bin/g++-4.9 && make install
#########

######### libpll-1.0.2
WORKDIR /opt/
RUN wget http://www.libpll.org/Downloads/libpll-1.0.2-sse3-64.tar.gz && \
    tar xzf libpll-1.0.2-sse3-64.tar.gz && rm libpll-1.0.2-sse3-64.tar.gz
WORKDIR /opt/libpll-1.0.2-sse3-64/
RUN cp libpll-sse3* /usr/local/lib && cp -r include/pll /usr/local/include/

RUN git clone https://github.com/Boussau/PHYLDOG /usr/local/phyldog && \
    mkdir /usr/local/phyldog/build
WORKDIR /usr/local/phyldog/build

ENV LD_LIBRARY_PATH /usr/local/lib/
RUN cmake ..  &&  make install
#########

#### build a working directory
RUN mkdir -p /data
WORKDIR /data

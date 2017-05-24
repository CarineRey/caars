FROM carinerey/apytram
MAINTAINER Carine Rey carine.rey@ens-lyon.fr

RUN apt-get clean && apt-get update

#### install ete2
RUN apt-get install --no-install-recommends -qy unzip \
                        ruby \
#to install ete2
                        python-setuptools \
                        python-qt4 \
                        python-scipy \
                        python-mysqldb \
                        python-lxml \
                        python-pip \
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
                        g++-4.9 \
                        libboost-all-dev \
                        ssh \ 
#mpi needs something to communicate (rsh, ssh other)
# to install fasttree
                        fasttree


RUN pip install ete2

#### install TransDecoder
WORKDIR /opt
RUN wget https://github.com/TransDecoder/TransDecoder/archive/v3.0.0.zip
RUN unzip v3.0.0.zip  && rm v3.0.0.zip
WORKDIR /opt/TransDecoder-3.0.0/
RUN make
ENV PATH /opt/TransDecoder-3.0.0/:$PATH


### Install muscle
WORKDIR /usr/local/bin
RUN wget http://www.drive5.com/muscle/downloads3.8.31/muscle3.8.31_i86linux64.tar.gz && \
    tar xvzf muscle3.8.31_i86linux64.tar.gz && rm muscle3.8.31_i86linux64.tar.gz && mv muscle3.8.31_i86linux64 muscle


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

RUN gpg --keyserver ha.pool.sks-keyservers.net --recv-keys B42F6819007F00F88E364FD4036A9C25BF357DD4
RUN curl -o /usr/local/bin/gosu -SL "https://github.com/tianon/gosu/releases/download/1.4/gosu-$(dpkg --print-architecture)" \
    && curl -o /usr/local/bin/gosu.asc -SL "https://github.com/tianon/gosu/releases/download/1.4/gosu-$(dpkg --print-architecture).asc" \
    && gpg --verify /usr/local/bin/gosu.asc \
    && rm /usr/local/bin/gosu.asc \
    && chmod +x /usr/local/bin/gosu




### install opam
WORKDIR /home/user
ENV HOME /home/user
RUN wget https://raw.github.com/ocaml/opam/master/shell/opam_installer.sh -O - | sh -s /usr/local/bin 4.03.0
ENV PATH /usr/local/bin/.opam/4.03.0/bin:$PATH
RUN /usr/local/bin/opam init --comp 4.03.0
RUN echo ". /home/user/.opam/opam-init/init.sh > /dev/null 2> /dev/null || true" >> ~/.profile
#RUN echo "export PATH=/usr/local/bin/:$PATH" >> ~/.profile

ENV PATH /usr/local/bin/:/home/user/.opam/4.03.0/bin/:$PATH
ENV CAML_LD_LIBRARY_PATH /home/user/.opam/4.03.0/lib/stublibs
ENV OCAML_TOPLEVEL_PATH /home/user/.opam/4.03.0/lib/toplevel
ENV MANPATH /home/user/.opam/4.03.0/man
ENV PERL5LIB /home/user/.opam/4.03.0/lib/perl5

##### install bistro
RUN opam pin add bistro --dev-repo

### install amalgam
WORKDIR /opt/amalgam/
RUN git clone https://github.com/carinerey/amalgam /opt/amalgam/ && make
RUN cp /opt/amalgam/etc/entrypoint.sh /usr/local/bin/entrypoint.sh
#ADD entrypoint.sh /usr/local/bin/entrypoint.sh
ENV PATH /opt/amalgam/utils/bin:/opt/amalgam:$PATH
ENV PYTHONPATH /opt/amalgam/utils/lib

##### build a working directory
WORKDIR /data

##### update amalgam
WORKDIR /opt/amalgam/
RUN git fetch origin && git checkout test_aligner
RUN make
RUN cp /opt/amalgam/etc/entrypoint.sh /usr/local/bin/entrypoint.sh
#ADD entrypoint.sh /usr/local/bin/entrypoint.sh

RUN apt-get clean && apt-get update && apt-get install --no-install-recommends -qy sudo

CMD ["amalgam_app.byte", "-h"]
ENTRYPOINT ["bash", "/usr/local/bin/entrypoint.sh"]

FROM carinerey/apytram:v1.2-alpha
MAINTAINER Carine Rey carine.rey@ens-lyon.org

RUN apt-get clean && apt-get update

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
RUN wget https://raw.github.com/ocaml/opam/master/shell/opam_installer.sh -O - | sh -s /usr/local/bin 4.03.0
ENV PATH /usr/local/bin/.opam/4.03.0/bin:$PATH

RUN /usr/local/bin/opam init --comp 4.03.0
RUN echo ". $HOME/.opam/opam-init/init.sh > /dev/null 2> /dev/null || true" >> ~/.profile

ENV PATH $HOME/.opam/4.03.0/bin/:$PATH
ENV CAML_LD_LIBRARY_PATH $HOME/.opam/4.03.0/lib/stublibs
ENV OCAML_TOPLEVEL_PATH $HOME/.opam/4.03.0/lib/toplevel
ENV MANPATH $HOME/.opam/4.03.0/man:$MANPATH
ENV PERL5LIB $HOME/.opam/4.03.0/lib/perl5:$PERL5LIB

##### install bistro
#RUN opam install -y bistro=0.3.0
RUN opam pin add bistro --dev-repo

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

##### install entrypoint
ADD entrypoint.sh /usr/local/bin/entrypoint.sh

##### install local dependencies
WORKDIR /opt/caars_py_lib
ENV PYTHONPATH  /opt/caars_py_lib:$PYTHONPATH
COPY utils/lib/*.py /opt/caars_py_lib/
COPY utils/bin/* /usr/bin/


#### install treerecs
WORKDIR /opt/treerecs
RUN git clone https://gitlab.inria.fr/Phylophile/Treerecs.git .
RUN cmake -DCMAKE_BUILD_TYPE=MinSizeRel .
RUN make Treerecs
ENV PATH /opt/treerecs/bin/:$PATH


##### build a working directory
WORKDIR /data
ENTRYPOINT ["bash", "/usr/local/bin/entrypoint.sh"]

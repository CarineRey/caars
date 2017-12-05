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
RUN git clone --branch rel220  https://github.com/BioPP/bpp-core && \
    git clone --branch rel220  https://github.com/BioPP/bpp-seq && \
    git clone --branch rel220  https://github.com/BioPP/bpp-phyl
RUN cd  $bpp_dir/sources_bpp/bpp-core &&\
    cmake . &&\
    (make -j 4 || make) &&\
    make install &&\
    cd $bpp_dir/sources_bpp/bpp-seq &&\
    cmake  . &&\
    (make -j 4 || make) &&\
    make install &&\
    cd  $bpp_dir/sources_bpp/bpp-phyl &&\
    cmake - . &&\
    (make -j 4 || make) &&\
    make install &&\
    rm -r $bpp_dir/sources_bpp

##### install phylomerge
RUN git clone https://github.com/boussau/phylomerge/ /opt/phylomerge
WORKDIR /opt/phylomerge
RUN bash ./compile
ENV PATH /opt/phylomerge:$PATH

#### install profileNJ
RUN pip install --upgrade pip
RUN pip install ete3==3.0.0b35
RUN pip install profileNJ

##### install phyldog

# ######### libpll-1.0.11
# ##  WORKDIR /root/usr/bin
# ##  RUN wget http://www.libpll.org/Downloads/libpll-1.0.11.tar.gz && \
# ##      tar xzf libpll-1.0.11.tar.gz
# ##  
# ##  WORKDIR /root/usr/bin/libpll-1.0.11
# ##  ENV CXX /usr/bin/g++-4.9
# ##  RUN ./configure && make install
# ##  
# ##  WORKDIR /opt/phyldog/build
# ##  RUN cmake .. -DCMAKE_CXX_COMPILER=/usr/bin/g++-4.9 && make install
# #########
# 
# ######### libpll-1.0.2
# #WORKDIR /opt/
# #RUN wget http://www.libpll.org/Downloads/libpll-1.0.2-sse3-64.tar.gz && \
# #    tar xzf libpll-1.0.2-sse3-64.tar.gz && rm libpll-1.0.2-sse3-64.tar.gz
# #WORKDIR /opt/libpll-1.0.2-sse3-64/
# #RUN cp libpll-sse3* /usr/local/lib && cp -r include/pll /usr/local/include/


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
RUN opam install -y bistro=0.3.0

######### libpll (new depot github)
#WORKDIR /opt/
#RUN apt-get clean && apt-get update && apt-get install --no-install-recommends -qy flex bison autotools-dev autoconf libtool automake
#RUN git clone --branch 0.3.2 https://github.com/xflouris/libpll.git && \
#    cd libpll  && ./autogen.sh &&  ./configure && make && make install  
RUN opam pin add bistro --dev-repo

######### libpll-1.0.2
COPY libpll-1.0.2-sse3-64.tar.gz /opt/libpll-1.0.2-sse3-64.tar.gz
WORKDIR /opt/
RUN tar xzf libpll-1.0.2-sse3-64.tar.gz && rm libpll-1.0.2-sse3-64.tar.gz
WORKDIR /opt/libpll-1.0.2-sse3-64/
RUN cp libpll-sse3* /usr/local/lib && cp -r include/pll /usr/local/include/


######### phyldog
# for libpll
RUN apt-get clean && apt-get update && apt-get install --no-install-recommends -qy flex bison autotools-dev autoconf libtool automake

WORKDIR  /opt/phyldog/

RUN git clone --branch faster_reconciliation  https://github.com/BenoitMorel/PHYLDOG.git /opt/phyldog/ && \
    git checkout 938e6b94bb659db5ffcdfefcf97420b7decb8c4b

RUN sed -e "s|git@github.com:ddarriba/pll-modules.git| --branch dev https://github.com/ddarriba/pll-modules.git|" -i ./build_libpll2.sh && ./build_libpll2.sh

RUN mv /opt/phyldog/libpll2_dependencies/include/pll /opt/phyldog/libpll2_dependencies/include/pllmodules
RUN mkdir -p build && cd /opt/phyldog/build && cmake  -DCMAKE_LIBRARY_PATH=/opt/phyldog/libpll2_dependencies/libs -DCMAKE_INCLUDE_PATH=/opt/phyldog/libpll2_dependencies/include ..  &&  make -j 4 install
ENV PATH /opt/phyldog/build/bin/:$PATH
#########

## update apytram
WORKDIR /opt/apytram
RUN git pull && git checkout ab9e8eb

##### install entrypoint
ADD entrypoint.sh /usr/local/bin/entrypoint.sh


##### install local dependencies
WORKDIR /opt/caars_py_lib
COPY utils/bin/* /usr/bin/
COPY utils/lib/*.py /opt/caars_py_lib/
ENV PYTHONPATH  /opt/caars_py_lib:$PYTHONPATH

##### build a working directory
WORKDIR /data
ENTRYPOINT ["bash", "/usr/local/bin/entrypoint.sh"]

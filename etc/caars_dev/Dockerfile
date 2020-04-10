FROM carinerey/caars_env:generax_20200331
MAINTAINER Carine Rey carine.rey@ens-lyon.org

### install caars master
WORKDIR /opt/caars/
RUN git clone https://github.com/carinerey/caars /opt/caars/
ENV PATH /opt/caars/utils/bin:/opt/caars:$PATH
ENV PYTHONPATH /opt/caars/utils/lib

##### update caars on the dev branch
ARG BRANCH_DEV
ENV BRANCH_DEV=${BRANCH_DEV}
RUN echo BRANCH_DEV $BRANCH_DEV && git fetch origin && git checkout $BRANCH_DEV
RUN make && (make test || make test) && make clean_test

##### build a working directory
WORKDIR $HOME

CMD ["caars", "-h"]

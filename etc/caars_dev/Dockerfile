ARG ENV_TAG
FROM carinerey/caars_env:${ENV_TAG}
MAINTAINER Carine Rey carine.rey@ens-lyon.org

### install caars master
WORKDIR /opt/caars/
RUN git clone https://github.com/carinerey/caars /opt/caars/
ENV PATH /opt/caars/utils/bin:/opt/caars:$PATH
ENV PYTHONPATH /opt/caars/utils/lib

##### update caars on the dev branch
ARG BRANCH
ENV BRANCH=${BRANCH}
RUN echo BRANCH $BRANCH && git fetch origin && git checkout $BRANCH
RUN make && (make test || make test) && make clean_test

##### build a working directory
WORKDIR $HOME

CMD ["caars", "-h"]

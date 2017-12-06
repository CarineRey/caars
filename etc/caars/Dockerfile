FROM carinerey/caars_env:0.2.0-bistro3
MAINTAINER Carine Rey carine.rey@ens-lyon.org

### install caars master
WORKDIR /opt/caars/
RUN git clone https://github.com/carinerey/caars /opt/caars/
ENV PATH /opt/caars/utils/bin:/opt/caars:$PATH
ENV PYTHONPATH /opt/caars/utils/lib

RUN make && (make test || make test) && make clean_test
##### build a working directory
WORKDIR $HOME

CMD ["caars", "-h"]


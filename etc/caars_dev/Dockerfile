FROM carinerey/apytramFROM carinerey/caars_env
MAINTAINER Carine Rey carine.rey@ens-lyon.org

### install caars master
WORKDIR /opt/caars/
RUN git clone https://github.com/carinerey/caars /opt/caars/
ENV PATH /opt/caars/utils/bin:/opt/caars:$PATH
ENV PYTHONPATH /opt/caars/utils/lib

##### update caars on the dev branch
WORKDIR /opt/caars/
RUN git fetch origin && git checkout update-wrt-bistro
RUN make clean && make

##### build a working directory
WORKDIR /data

CMD ["caars", "-h"]
ENTRYPOINT ["bash", "/usr/local/bin/entrypoint.sh"]

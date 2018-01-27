FROM ubuntu:16.04

WORKDIR /data/

COPY docker-build.sh .
COPY lumpy-smoother .

RUN bash docker-build.sh

ENTRYPOINT ["lumpy-smoother"]

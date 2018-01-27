FROM ubuntu:17.04

WORKDIR /data/

COPY docker-build.sh .

RUN bash docker-build.sh

ENTRYPOINT ["lumpy-smoother"]

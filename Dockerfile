FROM ubuntu:16.04


COPY docker-build.sh .
COPY lumpy-smoother .

RUN bash docker-build.sh

WORKDIR /work/

#ENTRYPOINT ["lumpy-smoother"]

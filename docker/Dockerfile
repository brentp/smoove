FROM ubuntu:16.04

COPY docker-build.sh .
COPY ./smoove .

RUN bash -x docker-build.sh

WORKDIR /work/

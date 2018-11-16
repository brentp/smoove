FROM ubuntu:16.04

COPY docker-build.sh .
COPY ./smoove .

RUN touch asdf.xxx
RUN bash docker-build.sh
RUN rm -f asdf.xxx

WORKDIR /work/

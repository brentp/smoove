FROM ubuntu:16.04

RUN apt update && \
    apt-get install -qy software-properties-common &&  \
    apt-add-repository ppa:deadsnakes/ppa && apt update

RUN apt-get -qy install python3.7 python3.7-dev curl && \
    curl https://bootstrap.pypa.io/get-pip.py -o get-pip.py && \
    python3.7 get-pip.py

COPY docker-build.sh .
COPY ./smoove .

RUN touch asdf.xxx
RUN bash docker-build.sh
RUN rm -f asdf.xxx

WORKDIR /work/

FROM ubuntu:16.04
ENV PATH=/root/miniconda2/bin/:$PATH

RUN apt-get update && \
    apt-get -qy install wget curl git \
    zlib1g-dev \
    make build-essential cmake libncurses-dev ncurses-dev g++ gcc \
    nfs-common \
    pigz bedtools gawk curl fuse wget git mdadm time \
    libbz2-dev lzma-dev liblzma-dev \
    syslog-ng libssl-dev libtool autoconf automake \
    libcurl4-openssl-dev libffi-dev libblas-dev liblapack-dev libatlas-base-dev

RUN wget -q https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh && \
  bash Miniconda2-latest-Linux-x86_64.sh -fb

COPY ./docker-build.sh .

RUN touch asdf.xxx
RUN bash docker-build.sh
RUN rm -f asdf.xxx

COPY ./smoove /usr/local/bin/smoove

WORKDIR /work/

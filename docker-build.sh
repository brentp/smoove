#!/bin/bash

set -euo pipefail

# used by Dockerfile
apt-get update
apt-get -qy install \
    zlib1g-dev \
    make build-essential cmake libncurses-dev ncurses-dev g++ gcc \
    python python-dev python-pip nfs-common \
    pigz bedtools gawk curl fuse wget git mdadm time \
    libbz2-dev lzma-dev liblzma-dev \
    syslog-ng libssl-dev libtool autoconf automake \
    libcurl4-openssl-dev libffi-dev libblas-dev liblapack-dev libatlas-base-dev \

git clone --recursive https://github.com/samtools/htslib.git \
git clone --recursive https://github.com/samtools/samtools.git \
git clone --recursive https://github.com/samtools/bcftools.git \
cd htslib && git checkout 1.6 && autoheader && autoconf && ./configure --enable-libcurl \
cd .. && make -j4 -C htslib install \
cd samtools && git checkout 1.6 \
autoreconf && ./configure && make -j4 install \
cd .. && cp ./samtools/samtools /usr/local/bin/ \
pip install -U awscli cython slurmpy toolshed awscli-cwlogs pyvcf pyfaidx cyvcf2 \


git clone --single-branch --recursive --depth 1 https://github.com/arq5x/lumpy-sv
cd lumpy-sv
make
cp ./bin/* /usr/local/bin/


## CNVnator stuffs
basedir=$(pwd)
git clone --depth 1 http://root.cern.ch/git/root.git
mkdir root/ibuild
cd root/ibuild
cmake -D x11=OFF ../
make -j4 install
cd last
git clone --depth 1 https://github.com/abyzovlab/CNVnator
cd CNVnator
cp -r ../samtools/ .
cp -r ../htslib/ .
make -j4 HTSDIR=htslib/ LIBS="-llzma -lbz2 -lz -lcurl -lssl -lcrypto"
cp ./cnvnator /usr/local/bin


rm -rf /var/lib/apt/lists/*

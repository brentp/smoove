#!/bin/bash

set -euo pipefail
basedir=$(pwd)

if [[ -e ./smoove ]]; then
cp ./smoove /usr/local/bin
chmod +x /usr/local/bin/smoove
else
    # TODO
    echo "NotImplemented"
fi

# used by Dockerfile
apt-get update
apt-get -qy install \
    zlib1g-dev \
    make build-essential cmake libncurses-dev ncurses-dev g++ gcc \
    python python-dev python-pip nfs-common \
    pigz bedtools gawk curl fuse wget git mdadm time \
    libbz2-dev lzma-dev liblzma-dev \
    syslog-ng libssl-dev libtool autoconf automake \
    libcurl4-openssl-dev libffi-dev libblas-dev liblapack-dev libatlas-base-dev

git clone --depth 1 https://github.com/ebiggers/libdeflate.git 
cd libdeflate
make -j 2 CFLAGS='-fPIC -O3' libdeflate.a
cp libdeflate.a /usr/local/lib
cp libdeflate.h /usr/local/include
cd $basedir
rm -rf libdeflate

git clone --recursive https://github.com/samtools/htslib.git
git clone --recursive https://github.com/samtools/samtools.git
git clone --recursive https://github.com/samtools/bcftools.git
cd htslib && git checkout 1.9 && autoheader && autoconf && ./configure --enable-libcurl --with-libdeflate
cd .. && make -j4 CFLAGS="-fPIC -O3" -C htslib install
cd $basedir

cd bcftools && git checkout 1.9
autoreconf && ./configure
set +e
make bcftools "PLUGINS_ENABLED=no" #
#"CFLAGS=-g -Wall -O2 -pedantic -std=c99 -D_XOPEN_SOURCE=600"
set -e
cp ./bcftools /usr/local/bin
cd $basedir
rm -rf bcftools

cd samtools && git checkout 1.9
autoreconf && ./configure && make -j2 CFLAGS='-fPIC -O3' install
cd $basedir && cp ./samtools/samtools /usr/local/bin/

wget -qO /usr/bin/batchit https://github.com/base2genomics/batchit/releases/download/v0.4.2/batchit
chmod +x /usr/bin/batchit

export HTSLIB_LIBRARY_DIR=/usr/local/lib
export HTSLIB_INCLUDE_DIR=/usr/local/include
pip install numpy pysam awscli cython toolshed awscli-cwlogs pyvcf pyfaidx cyvcf2 pip


cd $basedir
git clone -b no-big-ci https://github.com/brentp/svtyper
cd svtyper && python setup.py install
cd $basedir
rm -rf svtyper

cd $basedir
git clone https://github.com/hall-lab/svtools
cd svtools && python setup.py install
cd $basedir
rm -rf svtools

wget -qO /usr/local/bin/mosdepth https://github.com/brentp/mosdepth/releases/download/v0.2.3/mosdepth
chmod +x /usr/local/bin/mosdepth

wget -qO /usr/local/bin/duphold https://github.com/brentp/duphold/releases/download/v0.0.7/duphold
chmod +x /usr/local/bin/duphold

wget -qO /usr/bin/gsort https://github.com/brentp/gsort/releases/download/v0.0.6/gsort_linux_amd64
chmod +x /usr/bin/gsort

wget -qO /usr/bin/gargs https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux
chmod +x /usr/bin/gargs

git clone --single-branch --recursive --depth 1 https://github.com/arq5x/lumpy-sv
cd lumpy-sv
make -j 3
cp ./bin/* /usr/local/bin/

cd $basedir

rm -rf lumpy-sv

ldconfig

rm -rf /var/lib/apt/lists/*

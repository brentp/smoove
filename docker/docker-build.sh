#!/bin/bash

set -euo pipefail
basedir=$(pwd)

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

## CNVNATOR
export CPLUS_INCLUDE_PATH=/usr/include/root/

cd $basedir
git clone --depth 1 -b stdin https://github.com/brentp/CNVnator
cd CNVnator
ln -s $basedir/samtools/ .
ln -s $basedir/htslib/ .
make -j4 HTSDIR=htslib/ LIBS="-llzma -lbz2 -lz -lcurl -lssl -lcrypto -ldeflate"

cp ./cnvnator /usr/local/bin
cd $basedir
rm -rf CNVnator

wget -qO /usr/bin/batchit https://github.com/base2genomics/batchit/releases/download/v0.4.2/batchit
chmod +x /usr/bin/batchit

export HTSLIB_LIBRARY_DIR=/usr/local/lib
export HTSLIB_INCLUDE_DIR=/usr/local/include
export LD_LIBRARY_PATH=/usr/local/lib

wget -qO /usr/local/bin/mosdepth https://github.com/brentp/mosdepth/releases/download/v0.2.6/mosdepth
chmod +x /usr/local/bin/mosdepth

wget -qO /usr/local/bin/duphold https://github.com/brentp/duphold/releases/download/v0.2.1/duphold
chmod +x /usr/local/bin/duphold

wget -qO /usr/bin/gsort https://github.com/brentp/gsort/releases/download/v0.0.6/gsort_linux_amd64
chmod +x /usr/bin/gsort

wget -qO /usr/bin/bpbio https://github.com/brentp/bpbio/releases/download/v0.1.7/bpbio_static
chmod +x /usr/bin/bpbio

wget -qO /usr/bin/gargs https://github.com/brentp/gargs/releases/download/v0.3.9/gargs_linux
chmod +x /usr/bin/gargs

wget -qO /usr/bin/goleft https://github.com/brentp/goleft/releases/download/v0.2.1/goleft_linux64
chmod +x /usr/bin/goleft

wget -qO /usr/bin/smoove https://github.com/brentp/smoove/releases/download/v0.2.4/smoove
chmod +x /usr/bin/smoove

git clone --single-branch --recursive --depth 1 https://github.com/arq5x/lumpy-sv
cd lumpy-sv
make -j 3
cp ./bin/* /usr/local/bin/

cd $basedir

rm -rf lumpy-sv
rm -rf /var/lib/apt/lists/*

echo "done"

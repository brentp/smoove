#!/bin/bash

cp ./smoove $(which smoove)

rm -f to
mkdir to
REF=tests/subset.fa.gz
export TMPDIR=`pwd` && smoove call -o to --processes 2 -x --genotype --fasta $REF --name NA24385-svs-cram --excludechroms '~^GL,~^HLA,~_random,~^chrUn,~alt,~decoy' tests/NA24385_chr1.cram

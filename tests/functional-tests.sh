#!/bin/bash

set -eu
cp ./smoove $(which smoove)

rm -rf to
mkdir to
REF=tests/subset.fa.gz
<<DONE
export TMPDIR=`pwd` && smoove call -o to --processes 2 -x --genotype --fasta $REF --name NA24385-svs-cram --excludechroms '~^GL,~^HLA,~_random,~^chrUn,~alt,~decoy' tests/NA24385_chr1.cram

smoove call -o to --processes 2 --fasta $REF --name NA24385-svs-cram --excludechroms '~^GL,~^HLA,~_random,~^chrUn,~alt,~decoy' tests/NA24385_chr1.cram

smoove merge --fasta $REF -o to --name NA24385-svs-cram to/NA24385-svs-cram-smoove.vcf.gz
smoove genotype -x --fasta $REF -o to --name NA24385-svs-cram-g --vcf to/NA24385-svs-cram.sites.vcf.gz tests/NA24385_chr1.cram 
DONE

test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest

. ssshtest


run check_smoove_annotate go run cmd/smoove/smoove.go annotate --gff tests/test-anno.gff tests/test-anno.vcf
assert_exit_code 0
assert_equal 1 $(tail -1 $STDOUT_FILE | grep -c "exon:1:749")

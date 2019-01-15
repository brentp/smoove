# smoove  

[![Build Status](https://travis-ci.org/brentp/smoove.svg?branch=master)](https://travis-ci.org/brentp/smoove)

`smoove` simplifies and speeds calling and genotyping SVs for short reads. It also improves specificity by removing many
spurious alignment signals that are indicative of low-level noise and often contribute to spurious calls.

There is a blog-post describing `smoove` in more detail [here](https://brentp.github.io/post/smoove/)

It both supports small cohorts in a single command, and population-level calling with 4 total steps, 2
of which are parallel by sample.

It requires:

 + [lumpy and lumpy\_filter](https://github.com/arq5x/lumpy-sv)
 + [samtools](https://github.com/samtools/samtools): for CRAM support
 + [gsort](https://github.com/brentp/gsort): to sort final VCF
 + [bgzip+tabix](https://github.com/samtools/htslib): to compress and index final VCF

 And optionally (but all highly recommended):

 + [svtyper](https://github.com/hall-lab/svtyper): to genotypes SVs
 + [svtools](https://github.com/hall-lab/svtools): required for large cohorts
 + [mosdepth](https://github.com/brentp/mosdepth): remove high coverage regions.
 + [bcftools](https://github.com/samtools/bcftools): version 1.5 or higher for VCF indexing and filtering. 
 + [duphold](https://github.com/brentp/duphold): to annotate depth changes within events and at the break-points.

 Running `smoove` without any arguments will show which of these are found so they can be added to the PATH as needed.

`smoove` will:

1. parallelize calls to `lumpy_filter` to extract split and discordant reads required by lumpy
2. further filter `lumpy_filter` calls to remove high-coverage, spurious regions and user-specified chroms like 'hs37d5';
   it will also remove reads that we've found are likely spurious signals. 
   after this, it will remove singleton reads (where the mate was removed by one of the previous filters) from the discordant
   bams. This makes `lumpy` much faster and less memory-hungry.
3. calculate per-sample metrics for mean, standard deviation, and distribution of insert size as required by lumpy.
4. stream output of lumpy directly into multiple svtyper processes for parallel-by-region genotyping while lumpy is still running.
5. sort, compress, and index final VCF.

# installation

you can get `smoove` and all dependencies via (a large) docker image:

```
docker pull brentp/smoove
docker run -it brentp/smoove smoove -h
```

Or, you can download a `smoove` binary from here: https://github.com/brentp/smoove/releases
When run without any arguments, `smoove` will show you which of it's dependencies it can find
so you can adjust your $PATH and install accordingly.

# usage

## small cohorts (n < ~ 40)

for small cohorts it's possible to get a jointly-called, genotyped VCF in a **single command**.

```
smoove call -x --name my-cohort --exclude $bed --fasta $fasta -p $threads --genotype /path/to/*.bam
```
output will go to `./my-cohort-smoove.genotyped.vcf.gz`

the `--exclude $bed` is highly recommended as it can be used to ignore reads that overlap problematic regions.

A good set of regions for GRCh37 is [here](https://github.com/hall-lab/speedseq/blob/master/annotations/ceph18.b37.lumpy.exclude.2014-01-15.bed).

And for hg38 [here](https://github.com/hall-lab/speedseq/blob/master/annotations/exclude.cnvnator_100bp.GRCh38.20170403.bed)

## population calling

For population-level calling (large cohorts) the steps are:

1. For each sample, call genotypes:

```
smoove call --outdir results-smoove/ --exclude $bed --name $sample --fasta $fasta -p 1 --genotype /path/to/$sample.bam
```

For large cohorts, it's better to parallelize across samples rather than using a large $threads per sample. `smoove` can only
parallelize up to 2 or 3 threads on a single-sample and it's most efficient to use 1 thread.

output will go to `results-smoove/$sample-smoove.genotyped.vcf.gz``

2. Get the union of sites across all samples (this can parallelize this across as many CPUs or machines as needed):

```
# this will create ./merged.sites.vcf.gz
smoove merge --name merged -f $fasta --outdir ./ results-smoove/*.genotyped.vcf.gz
```

3. genotype all samples at those sites (this can parallelize this across as many CPUs or machines as needed) and run [duphold](https://github.com/brentp/duphold) to add depth annotations.

```
smoove genotype -d -x -p 1 --name $sample-joint --outdir results-genotped/ --fasta $fasta --vcf merged.sites.vcf.gz /path/to/$sample.$bam
```

4. paste all the single sample VCFs with the same number of variants to get a single, squared, joint-called file.

```
smoove paste --name $cohort results-genotyped/*.vcf.gz
```

5. (optional) annotate the variants with exons, UTRs that overlap from a GFF and annotate high-quality heterozygotes:

```
smoove annotate --gff Homo_sapiens.GRCh37.82.gff3.gz $cohort.smoove.square.vcf.gz | bgzip -c > $cohort.smoove.square.anno.vcf.gz
```

This adds a `SHQ` (Smoove Het Quality) tag to every sample format) a value of **4 is a high quality** call and the value of 1 is low quality. -1 is non-het.
It also adds a `MSHQ` for Mean SHQ to the INFO field which is the mean SHQ score across all heterozygous samples for that variant.

As a first pass, users can look for variants with MSHQ > 3. If you added [duphold](https://github.com/brentp/duphold) annotations, it's also
useful to check deletions with `DHFFC < 0.7` and duplications with `DHFFC > 1.25`.

# Troubleshooting

+ A panic with a message like ` Segmentation fault      (core dumped) | bcftools view -O z -c 1 -o` is likely to mean you have an old version of bcftools. 
  see #10

+ `smoove` will write to the system TMPDIR. For large cohorts, make sure to set this to something with a lot of space. e.g. `export TMPDIR=/path/to/big`

+ `smoove` requires recent version of `lumpy` and `lumpy_filter` so build those from source or get the most recent bioconda version.

# see also

[svtools](https://github.com/hall-lab/svtools)

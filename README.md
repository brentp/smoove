# smoove

`smoove` simplifies and speeds calling and genotyping SVs.

It requires:

 + [lumpy and lumpy\_filter](https://github.com/arq5x/lumpy-sv)

 And optionally:

 + [svtyper](https://github.com/hall-lab/svtyper): to genotypes SVs
 + [svtools](https://github.com/hall-lab/svtools): required for large cohorts
 + [samtools](https://github.com/samtools/samtools): for CRAM support
 + [gsort](https://github.com/brentp/gsort): to sort final VCF
 + [bgzip+tabix](https://github.com/samtools/htslib): to compress and index final VCF
 + [mosdepth](https://github.com/brentp/mosdepth)

 Running `smoove` without any arguments will show which of these are found so they can be added to the PATH as needed.

`smoove` will:

1. parallelize calls to `lumpy_filter` to extract split and discordant reads required by lumpy
2. further filter `lumpy_filter` calls to remove high-coverage, spurious regions and user-specified chroms like 'hs37d5';
   it will also remove reads that we've found are likely spurious signals. 
   after this, it will remove singleton reads (where the mate was removed by one of the previous filters) from the discordant
   bams. This makes `lumpy` much faster and less memory-hungry.
3. calculate per-sample metrics for mean, standard deviation, and distribution of insert size as required by lumpy.
4. correct the reference allele (lumpy always puts 'N')
5. stream output of lumpy directly into multiple svtyper processes for parallel-by-region genotyping while lumpy is still running.
6. sort, compress, and index final VCF.

# installation

I will release a binary shortly, meanwhile, you can get this and all dependencies via (a large) docker image:

```
docker pull brentp/smoove
docker run -it brentp/smoove -h
```

# usage

## small cohorts (n < ~ 40)

for small cohorts it's possible to get a jointly-called, genotyped VCF in a single command.

```
smoove call --name my-cohort --fasta $fasta -p $threads --genotype /path/to/*.bam > /dev/null
```

## large cohorts

For large cohorts the steps are:

1. For each sample, call genotypes:

```
smoove call --outdir results-smoove/ --name $sample --fasta $fasta -p $threads --genotype /path/to/$sample.bam > /dev/null
```

2. Get the union of sites across all samples:

```
# this will create ./merged.sites.vcf.gz
smoove merge --name merged -f $fasta --outdir ./ results-smoove/*.vcf.gz
```

3. genotype all samples at those sites.

```
smoove genotype -p 1 --name $sample-joint --outdir results-genotped/ --fasta $fasta --vcf merged.sites.vcf.gz /path/to/$sample.$bam
```

4. paste all the single sample VCFs with the same number of variants to get a single, squared, joint-called file.

```
smoove paste --name $cohort results-genotyped/*.vcf.gz
```

# TODO

+ [ ] annotate high-quality calls

# see also

[svtools](https://github.com/hall-lab/svtools)

# lumpy-smoother

`lumpy-smoother` simplifies and speeds running [lumpy](https://github.com/arq5x/lumpy-sv) and associated pre and post-processing steps.

It requires:

 + [lumpy and lumpy\_filter](https://github.com/arq5x/lumpy-sv)

 And optionally:

 + [cnvnator](https://github.com/abyzovlab/CNVnator): makes per-sample CNV calls that lumpy can use
 + [svtyper](https://github.com/hall-lab/svtyper): to genotypes SVs
 + [samtools](https://github.com/samtools/samtools): for CRAM support
 + [gsort](https://github.com/brentp/gsort): to sort final VCF
 + [bgzip+tabix](https://github.com/samtools/htslib): to compress and index final VCF

 Running `lumpy-smoother -h` will show which of these are found so they can be added to the PATH as needed.

`lumpy-smoother` will:

1. parallelize calls to `lumpy_filter` to extract split and discordant reads required by lumpy
2. further filter `lumpy_filter` calls to remove high-coverage, spurious regions and user-specified chroms like 'hs37d5';
   after this, it will remove singleton reads (where the mate was removed by one of the previous filters) from the discordant
   bams. This makes `lumpy` much faster and less memory-hungry.
3. parallelize calling cnvnator if it is on the $PATH, including splitting the reference as it requires.
   calls to `lumpy_filter` and `cnvnator` are in the same process-pool for maximum efficiency
4. calculate per-sample metrics for mean, standard deviation, and distribution of insert size as required by lumpy.
5. correct the reference allele (lumpy always puts 'N')
6. stream output of lumpy directly into multiple svtyper processes for parallel-by-region genotyping while lumpy is still running.
7. sort, compress, and index final VCF.

# installation

I will release a binary shortly, meanwhile, you can get this and all dependencies via (a large) docker image:

```
docker pull brentp/lumpy-smoother
docker run -it brentp/lumpy-smoother -h
```

# usage

run `lumpy-smoother -h` for full usage

```
lumpy-smoother \
        -n my-project \                        # arbitrary project name for file prefixes
        -f $reference \
        --outdir results \
        --processes 10 \                       # parallelize with this many processors.
        -C hs37d5,GL000192.1 \                 # comma-delimited list of chroms to exclude
        --exclude low-complexity-regions.bed \ # see: https://github.com/hall-lab/speedseq/tree/master/annotations 
        data/*.bam                             # any number of BAM or CRAM files
```

# TODO

+ [ ] annotate high-quality calls
+ [ ] (unlikely) isolate steps so that users can call, e.g.: 
    lumpy-smoother cnvs
    lumpy-smoother filter
    lumpy-smoother lumpy
    lumpy-smoother call

# limitations

this is limited to cohorts of ~ 30-40 or so.

# see also

[svtools](https://github.com/hall-lab/svtools)

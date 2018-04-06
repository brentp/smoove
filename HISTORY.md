v0.1.5
======
+ sensitivity increases for variants without split reads.
+ use a fixed length of 4 in merge (gives much better merging)
+ add new annotate command which takes a gff and annotates with gene names
+ increase max_reads argument to svtyper to reduce false negatives
+ remove QCFail and duplicate reads from split and disc.bams. This is now
  part of `smoove` but also fixed in `lumpy_filter` upstream.

v0.1.4
======
+ fix bug when using `-x/--removepr` with `--genotype` (thanks @brad for test-case)

v0.1.3
======
+ update command to match latest mosdepth

v0.1.2
======
+ bcftools and gsort are now required.
+ output never goes to STDOUT
+ output contigs to header
+ add svtools to docker
+ add option to remove PRPOS and PREND

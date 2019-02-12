v0.2.4
======
+ add ##reference=$fasta to vcf header. (#58)

v0.2.3
======
+ save disc and split counts and output a plot as QC
+ annotate: annotate with up and downstream 5KB. so any gene within 5kb will be annotated.
+ smoove merge now outputs a plot that shows the number of split and discordant reads for each sample.
+ smoove now filters about 25% more split reads by looking at how the read is consumed
  this results in a 5% reduction in false positive genome-in-a-bottle deletions > 1KB while leaving
  the true positives unchanged relative to version 0.2.2. (Thanks @davemcg for feedback on this feature).

v0.2.2
======
+ huge reduction in number of discordant reads sent to lumpy for some problematic samples.
  this is work by @ernfrid which reduces lumpy runtime for all cases and dramatically lowers
  it for cases that were previously problematic due to large numbers of discordant reads.
  On our Genome in a Bottle test sample, this **drops the number of discordant reads
  from 1,434,810 to 332,384; a 4.3 fold reduction** while leaving the final smoove output unchanged
  except for 2 BNDs that were removed. This reduction in reads will be even more dramatic
  in problematic samples.
+ switch to use svtyper>=0.7.0 with `--max_ci_dist` parameter. To avoid this set: 
  `export SMOOVE_NO_MAX_CI=xx` # any value will work; however, it's higly recommended
  to use as-is.
+ cnvnator changes and filtering. (still requires github.com/brentp/CNVnator fork).
+ bugfix for case when no output directory was specified.
+ use mosdepth with --fast-mode **NOTE** this requires the latest version of mosdepth (0.2.4 or greater)

v0.2.1
======
+ fix bug in smoove duphold for samples > threads that resulted in stalling
+ smoove is now more discerning about reads that are soft-clipped on both ends as
  these could be due to inversions. if a read is not flipped relative to its mate
  and it is soft-clipped on both ends, it is still removed.
+ in some cases reads with a large NM (mismatches) were not filtered because they
  had an unexpected type (uint8). `smoove` now checks more types and for some bams
  will remove more records which improves specificity.
+ smoove will now be more conservative with `NM` counting. an insertion of 7 bases
  is counted as an NM of 7. smoove will now correct this to an NM of 1 so it counts
  the number of events rather than the number of bases of each event.
+ simplify rules for filtering and make them less strict.

v0.2.0
======
+ **add [duphold](https://github.com/brentp/duphold)** annotations. (requires duphold v0.0.4)

v0.1.11
=======
+ flip start and end for small percent of cases where that's a problem
+ expose min sample weight to allow adjusting required support for each variant

v0.1.10
======= 
+ fix type in annotate (s/Number/Integer/)

v0.1.9
======
+ better error message in merge.
+ use set -e in process calls so we don't exit with 0 even on failure

v0.1.8
======
+ don't use --sum in merge step. this reduces the size of CI when using `smoove merge`

v0.1.7
======
+ add a `SHQ`: `smoove het-quality` score to the FORMAT field and a `MSHQ`:`mean smoove het-quality` score
  to the INFO. Variants with a het-quality of 4 are quite good.

v0.1.6
======
+ sensitivity improvement by dropping reads with > 5 mismatches.

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

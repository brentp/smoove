v0.1.5
======
+ use a fixed length of 4 in merge (gives much better merging)
+ add new annotate command which takes a gff

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

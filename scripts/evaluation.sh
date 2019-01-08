set -euo pipefail

truvari=~/src/truvari/truvari.py
fasta=/data/human/g1k_v37_decoy.fa
cram=/data/human/hg002.cram

#ncftp -R ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/HG002_HiSeq300x_fastq/140528_D00360_0018_AH8VC6ADXX/ .

#for f in $(find . -type f -name "*_R1_*.fastq.gz" | sort); do zcat $f; done | bgzip -@ 14 -c > hg002_R1.fastq.gz &
#for f in $(find . -type f -name "*_R2_*.fastq.gz" | sort); do zcat $f; done | bgzip -@ 14 -c > hg002_R2.fastq.gz &
# ftp://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/

#wait
<<ALIGN
~/Projects/src/bwa/bwa mem -R '@RG\tID:HG002\tSM:HG002\tPL:ILLUMINA\tPU:HG002\tLB:HG002' -t 32 $fasta hg002_R1.fastq.gz hg002_R2.fastq.gz \
	| samblaster \
	| samtools sort -@ 4 -m 15G --output-fmt CRAM --reference $fasta -o hg002.cram

bcftools view -i SVTYPE="DEL" -O z -o HG002_SVs_Tier1_v0.6.DEL.vcf.gz HG002_SVs_Tier1_v0.6.vcf.gz
tabix HG002_SVs_Tier1_v0.6.DEL.vcf.gz
bcftools view -i 'REPTYPE="DUP"' -O z -o HG002_SVs_Tier1_v0.6.DUP.vcf.gz HG002_SVs_Tier1_v0.6.vcf.gz
tabix HG002_SVs_Tier1_v0.6.DUP.vcf.gz


samtools index hg002.cram
ALIGN

bcftools view -i 'SVTYPE="DEL"' -O z -o HG002_SVs_Tier1_v0.6.DEL.vcf.gz HG002_SVs_Tier1_v0.6.vcf.gz
tabix HG002_SVs_Tier1_v0.6.DEL.vcf.gz


#smoove_cmd=go run cmd/smoove/smoove.go
version=dev

rm -rf evaluation-$version
go run cmd/smoove/smoove.go call --genotype -o evaluation-$version/ -x -p 4 -f $fasta -n testing $cram

bcftools view -i 'SVTYPE="DEL"' evaluation-$version/testing-smoove.genotyped.vcf.gz -O z -o evaluation-$version/testing-smoove.genotyped.DEL.vcf.gz
tabix evaluation-$version/testing-smoove.genotyped.DEL.vcf.gz

ev=DUP
filt="> 1.3"

ev=DEL
filt="< 0.7"

odir=evaluation-$version/truvari

eval=$odir-eval-$ev
mkdir -p $eval/
rm -rf $eval/*

sizemax=15000000
sizemin=1000

poverlap=0.6

set -x
python $truvari --sizemax $sizemax -s $sizemin -S $((sizemin - 30)) -b HG002_SVs_Tier1_v0.6.$ev.vcf.gz -c evaluation-$version/testing-smoove.genotyped.$ev.vcf.gz -o $eval/smoove/ --passonly --pctsim=0  -r 20 --giabreport -f $fasta --no-ref --includebed HG002_SVs_Tier1_v0.6.bed -O $poverlap
python scripts/totable.py $eval/smoove/summary.txt

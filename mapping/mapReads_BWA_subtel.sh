#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 24:00:00
#SBATCH -J bwa_mapping

# bwa index ref.fa

module load bioinfo-tools
module load bwa/0.7.17

REF=$(realpath $1)
READS=$(realpath $2)
OUT=$(realpath $3)

# Script for mapping single-end reads to a reference genome using BWA
# Alternatively, the script mapped the reads to the hybrid hg38 + humanSTF500 for the analysis of subtelomeric reads (deprecated as the refined T2T human genome assembly is now available)

for READ in $(ls ${READS}*fq.gz)
do
	bwa aln ${REF} ${READS} > ${OUT}${READS/.fq.gz/_aln_sa.sai}
# The following allows to filter out reads with more than n multimapping positions, from the Stong et al. article it is proposed to allow multi mapped reads with up to 101 positions.
#bwa samse -n 101 ${REF} ${OUT}aln_sa.sai reads.fq.gz > aln-se.sam

done
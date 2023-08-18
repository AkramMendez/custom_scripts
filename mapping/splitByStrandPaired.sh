#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 02:00:00
#SBATCH -J split_bam


module load bioinfo-tools
module load samtools

bam_folder=$1


for bam in $(ls ${bam_folder}/*.bam)
do
	#
# 1. alignments of the second in pair if they map to the forward strand
#
samtools view -b -f 128 -F 16 ${bam} > fwd1.bam
samtools index fwd1.bam

# 2. alignments of the first in pair if their mate maps to the forward strand
samtools view -b -f 64 -F 32 ${bam} > fwd2.bam
samtools index fwd2.bam

#
# Combine alignments that originate on the forward strand.
#
samtools merge -f fwd.bam fwd1.bam fwd2.bam
samtools index fwd.bam

# 1. alignments of the second in pair if it maps to the reverse strand
#
samtools view -b -f 144 ${bam} > rev1.bam
samtools index rev1.bam

# 2. alignments of the first in pair if their mates map to the reverse strand
samtools view -b -f 96 ${bam} > rev2.bam
samtools index rev2.bam

#
# Combine alignments that originate on the reverse strand.
#
samtools merge -f rev.bam rev1.bam rev2.bam
samtools index rev.bam

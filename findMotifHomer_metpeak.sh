#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 12:00:00
#SBATCH -J computeMatrix

module load bioinfo-tools
module load HOMER/4.11

peaks=$(realpath $1) # Path to the directory with peak files (BED or narrowPeak format)
ref=$(realpath $2) # Path to the reference genome FASTA file.

# Script to analyze the sequences for motifs in regions containting m6A peaks (previously obatined with MetPeak or MACS2)

for peaks in $(ls ${peaks}/*.{bed}) 
do 
	sample=$(basename ${peaks})
	echo "${sample}" 
	findMotifsGenome.pl ${peaks} ${ref} homerMotifs_${sample%.bed} -size 75 -len 6 -basic -rna -p 8
done


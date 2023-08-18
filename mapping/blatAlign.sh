#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 08:00:00
#SBATCH -J blat_aln


hg38="GRCh38.primary_assembly.genome.fa"

subtel="humanSTF500.fasta"

echo "Aligning Subtelomere sequences to GRCh38 genome:"
blat ${hg38} ${subtel} GRCh38.subtel.blat.psl

echo "Alignment finished."

#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 00:30:00
#SBATCH -J igvtools


module load bioinfo-tools
module load IGVtools

gtf=$1

echo "Sorting GTF file:"

igvtools sort ${gtf} ${gtf%.gtf}.sorted.gtf

echo "Indexing sorted GTF file:"

igvtools index ${gtf%.gtf}.sorted.gtf

echo "Done."
#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 06:00:00
#SBATCH -J trimming_reads

module load bioinfo-tools
module load TrimGalore

rawreads=$(realpath $1)
outdir=$(realpath $2)

for i in $(ls ${rawreads}/*.fastq.gz);do trim_galore ${i} --cores  -o ${rawreads} ;done

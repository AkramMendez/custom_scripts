#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 10
#SBATCH -t 02:00:00
#SBATCH -J fastQC


module load bioinfo-tools
module load FastQC
module load MultiQC

reads_folder=$(realpath $1)
out=$(realpath $2)

echo "QC analysis (fastQC):"

fastqc ${reads_folder}/*.fastq.gz --outdir ${out} --threads ${SLURM_NTASKS}

echo "Merging fastQC reports (MultiQC):"

multiqc ${out} --outdir ${out}
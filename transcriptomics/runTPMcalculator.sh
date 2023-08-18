#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 10
#SBATCH -t 08:00:00
#SBATCH -J TPMcalc

conda activate xengsort
module load bioinfo-tools

inputdir=$(realpath $1)
gtf=$(realpath $2)

#This script calculates TPM values from paired-end BAM files located in a single directory using NCBI's TPMcalculator.
# For calculating single-end BAM files check the documentation.
#Inputdir = BAM files directory
# For a description of the resulting output files check: https://github.com/ncbi/TPMCalculator/wiki/Output-files

echo "Running TPMCalculator for BAM files in ${inputdir}"

TPMCalculator -g ${gtf} -k gene_name -t transcript_name -p -d ${inputdir}

echo "Done."


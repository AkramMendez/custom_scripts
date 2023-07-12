#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 6
#SBATCH -t 03:00:00
#SBATCH -J minionQC

module load bioinfo-tools
module load R_packages

minionQC="/domus/h1/amendez/private/minion_qc/MinIONQC.R"
inputdir=$(realpath $1) # Input directory path containg the 'sequencing_summary.txt' file
outdir=$(realpath $2) # Output directory path

# This script generates MinionQC (https://github.com/roblanf/minion_qc) reports for the quality control of Oxford Nanopore dRNA data
# , the reports are stored into a desired output directory
echo "Running minionQC"

Rscript ${minionQC} -i ${inputdir}/sequencing_summary.txt -o ${outdir} -p ${SLURM_NTASKS}

echo "Done."
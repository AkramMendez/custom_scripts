#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 18:00:00
#SBATCH -J metdiff

module load R_packages/4.0

# Script for implementing the m6A analysis pipeline with MetDiff (https://github.com/compgenomics/MeTDiff/) in SLURM
Rscript --vanilla --verbose metdiffAnalysis.R 2> logfile.runmetdiff.log

echo "Done."
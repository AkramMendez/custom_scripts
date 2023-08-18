#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 48:00:00
#SBATCH -J metdiff

module load R_packages/4.0

Rscript --vanilla --verbose metdiffAnalysis_infectedCells_samples_all.R 2> logfile.runmetdiff_infectedCells_samples_all.log

echo "Done."
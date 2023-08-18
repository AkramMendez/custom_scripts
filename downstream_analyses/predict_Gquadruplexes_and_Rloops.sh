#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 12:00:00
#SBATCH -J G4Find

fasta=$(realpath $1) # Reference genome fasta file
outdir=$(realpath $2) # Output directory

G4hunter="/domus/h1/amendez/private/G4-hunter/G4Hunter.py"
RloopFinder="/domus/h1/amendez/private/Rloop_finder.py"

# Script to predict R-loops and G-quadruplexes genome-wide using G4Hunter and Rloop finder tools

echo "Activating conda environment"
conda activate /crex/proj/nb_project/private/conda_envs/py2.7

cd ${outdir}

# Notes for default parameter definition from original author's publication: "Re-evaluation of G-quadruplex propensity with G4Hunter"
# "When analysing a genome-wide, the mean of the scored nucleic acid sequence is computed for a sliding window arbitrary set at 25 nt" from "Re-evaluation of G-quadruplex propensity with G4Hunter"
# "Based on the results of the analysis of the mitochondrial genome, a window size of 25 and a threshold of
# 1.5 results in precision above 90%, meaning that more that
# 90% of the sequences identified in the human genome with
# these settings should form G-quadruplexes"

echo "Predicting G-quadruplexes0"
python2.7 ${G4hunter} -i ${fasta} -o ${outdir} -w 25 -s 1.5

echo "Predicting R-loop forming sequences coordinates in BED format"
python2.7 ${RloopFinder} --model m1 -bed -i ${fasta} -o ${name}_RloopFinder_QmRLFS.bed --log --verbose

echo "Predicting R-loop forming sequences results in tabular format"
python2.7 ${RloopFinder} --model m1 -i ${fasta} -o ${name}_RloopFinder_QmRLFS_table.txt --log --verbose

echo "Done."



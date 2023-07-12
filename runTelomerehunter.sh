#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 12:00:00
#SBATCH -J telomerehunter

module load bioinfo-tools
module load TelomereHunter
module load samtools
module load R_packages/4.0

treatedbam=$(realpath $1)
controlbam=$(realpath $2)
samplename=$3
cytobands=$(realpath $4)
outdir=$(realpath $5)

# Script for analyzing the telomere content from matched 'treatment' and 'control' samples using TelomereHunter
# https://www.dkfz.de/en/applied-bioinformatics/telomerehunter/telomerehunter.html

echo "Running TelomereHunter"

telomerehunter -ibt ${treatedbam} \
-ibc ${controlbam} \
-p ${samplename} \
-b ${cytobands} \
--removeDuplicates \
--parallel \
--plotRevCompl --plotChr --plotFractions --plotTelContent --plotGC --plotRepeatFreq --plotTVR --plotSingleton --plotRevCompl
-o ${outdir}

echo "Done."
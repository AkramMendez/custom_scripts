#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 12:00:00
#SBATCH -J computeMatrix


module load bioinfo-tools
module load BEDTools

bigwig=$1
sample_name=$2

bed="/crex/proj/nb_storage/private/m6a_project/m6a_analysis_lab/csh/hisat2_alns/bamCoverageFiles/profiles.mettl3.csh/hg38.ensGene.genPred.bed"

computeMatrix scale-regions -S ${bigwig} -R ${bed} --beforeRegionStartLength 3000 --regionBodyLength 5000 --afterRegionStartLength 3000 --binSize 5 --skipZeros -o ${sample_name}.matrix.mat.gz -p ${SLURM_NTASKS} --outFileNameMatrix ${sample_name}.matrix.tab --verbose

mat=${sample_name}.matrix.mat.gz

echo "Plotting profile:"

plotProfile -m ${mat} -out ${sample_name}.profile.png 

echo "Done."
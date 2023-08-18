#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 12:00:00
#SBATCH -J computeMatrix


module load bioinfo-tools
module load deepTools

bigwigs_folder=$1
n=${SLURM_NTASKS}
ref="mm10.ensGene.bed"

cd ${TMPDIR}

echo "computing matrix for file(s): ${bigwigs}"

computeMatrix scale-regions -S $(ls ${bigwigs_folder}/*.bw) -R ${ref} \
--beforeRegionStartLength 3000 \
--regionBodyLength 5000 \
--afterRegionStartLength 3000 \
--skipZeros -o scaleRegions.log2ratio.matrix.mat.gz \
--outFileNameMatrix scaleRegions.matrix \
--numberOfProcessors=${n} \
--samplesLabel $(ls ${bigwigs_folder}/*.bw | grep -oP "(.*)(?=.hisat2.*)") \
--verbose

echo "Done."
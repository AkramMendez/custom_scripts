#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 12:00:00
#SBATCH -J bigwig

# Script to generate bigwig coverage tracks from BAM files for genomic visualization.

module load bioinfo-tools
module load deepTools/3.3.2

inputdir=$(realpath $1)
outdir=$(realpath $2)
label=$3

cd $SNIC_TMP

cp ${inputdir}/*.bam .
cp ${inputdir}/*.bai .

for bam in $(ls *.bam)
do 
	echo "Making Bigwig for file: ${bam}"
	bamCoverage -b ${bam} -p ${SLURM_NTASKS} --binSize 10 --outFileFormat bigwig -o ${bam%.bam}.${label}.bs10.bw
	
	echo "Saving bigwig to ${outdir}"
	cp ${bam%.bam}.${label}.bs10.bw ${outdir}

done

echo "Done."
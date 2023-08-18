#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 6
#SBATCH -t 03:00:00
#SBATCH -J trimming_reads

module load bioinfo-tools
module load TrimGalore/0.6.1

inputdir=$(realpath $1) # Input directory with raw reads in fastq.gz format
outdir=$(realpath $2) # Output directory

# Script for trimming paired-end raw reads for cleaning and adapter removal using Trimgalore

cd ${SNIC_TMP}

echo "Copying files to ${SNIC_TMP}"

cp ${inputdir}/*.fastq.gz ${SNIC_TMP}

echo "Trimming fastq files"
for i in $(ls *.fastq.gz | grep -oP "(.*)(?=\_\d+.fastq.gz)" | sort | uniq)
do 
	R1=$(echo "${i}_1.fastq.gz")
	R2=$(echo "${i}_2.fastq.gz") 
	
	echo "Processing file ${R1} and ${R2}"

	trim_galore ${R1} ${R2} --paired --length 20 --cores ${SLURM_NTASKS} -o ${outdir}
done 

echo "Done."

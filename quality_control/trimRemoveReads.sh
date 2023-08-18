#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 06:00:00
#SBATCH -J trimming_reads

module load bioinfo-tools
module load TrimGalore/0.6.1

rawreads=$(realpath $1)
outdir=$(realpath $2)

n=$(echo "$(($SLURM_NTASKS-1))")

for reads_file in $(ls ${rawreads}/*.fastq)
do 
	echo "Trimming reads: ${reads_file}"

	trim_galore ${reads_file} --cores ${n}  -o ${outdir}

	echo "Trimmed reads saved in ${outdir}"

	echo "Removing raw reads file"

	rm ${reads_file}

	echo "Done"

done
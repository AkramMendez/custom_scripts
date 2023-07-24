#!/bin/bash -l
#SBATCH -A naiss2023-22-153
#SBATCH -p core -n 10
#SBATCH -t 08:00:00
#SBATCH -J trimming_reads

module load bioinfo-tools
module load TrimGalore/0.6.1

inputdir=$(realpath $1) # Input directory with raw read FastQ files
outdir=$(realpath $2) # Output directory

cd ${SNIC_TMP}

echo "Copying files to ${SNIC_TMP}"

cp ${inputdir}/*.fastq.gz ${SNIC_TMP}

echo "Trimming fastq files"
for i in $(ls *.fastq.gz);do echo "Processing file ${i}" ;trim_galore ${i} --nextseq 20 --length 20 --cores ${SLURM_NTASKS} -o ${outdir} ;done 

echo "Done."
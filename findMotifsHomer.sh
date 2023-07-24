#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 03:00:00
#SBATCH -J findMotifs

# Script for finding motifs using HOMER from peak coordinates and a genome reference file.
module load bioinfo-tools
module load HOMER/4.11

inputdir=$(realpath $1) # Input directory with peak calling *.narrowPeak files
ref=$(realpath $2) # Reference genome FASTA file
name=$3 # Output name
outdir=$(realpath $4) # Output directory

cd ${SNIC_TMP}

cp ${inputdir}/*.narrowPeak .


findMotifs(){
	peaks=$1
	ref=$2
	name=$3
	outdir=$4

	findMotifsGenome.pl ${peaks} ${ref} ${outdir}/${name} -size 75 -len 6 -basic -rna -p ${SLURM_NTASKS}
}

for peaks in $(ls *.narrowPeak)
do 
	findMotifs ${peaks} ${ref} ${name} ${outdir}

done

echo "Done."
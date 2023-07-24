#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 12
#SBATCH -t 24:00:00
#SBATCH -J genMap

# Script to calculate genome-wide mappability scores suing GenMap https://github.com/cpockrandt/genmap. 
ref_fasta=$(realpath $1) # Reference genome FASTA file
index_dir_name=$2 # Directory name, genMap will create a directory after building the index
kmer_size=$3 #K-mer size can be estimated from average sequence length according to the specific sequencing data
mismatch=$4 #Allowed number of errors per k-mer
outdir=$(realpath $5)

module load bioinfo-tools

conda activate /proj/nb_project/private/conda_envs/genmap

echo "Building genMap index"

genmap index -F ${ref_fasta} -I ${index_dir_name}

echo "Computing mappability"

genmap -K ${kmer_size} \
-E ${mismatch} \
-I ${index_dir_name} \
-O ${outdir} \
-t -w -bg

echo "Done."
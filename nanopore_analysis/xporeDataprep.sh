#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 12
#SBATCH -t 48:00:00
#SBATCH -J xpdataprep

# Script to perform data preparation from Nanopore long-reads for m6A detection using xPore
# https://xpore.readthedocs.io/en/latest/quickstart.html
eventalign=$(realpath $1) # Path to eventaligned Nanopolish files
outdir=$(realpath $2) # Output directory
gtf=$(realpath $3) # Reference genome annotation GTF file
txm=$(realpath $4) # Reference transcriptome FASTA file

module load bioinfo-tools

echo "Activating xPore conda environment"

conda activate /proj/nb_project/private/conda_envs/nanopolish/

xpore dataprep --eventalign ${eventalign} \
--out_dir=${outdir} \
--gtf_or_gff=${gtf} \
--transcript_fasta=${txm} \
--chunk_size 10000 \
--genome \
--n_processes=${SLURM_NTASKS}

echo "Done."
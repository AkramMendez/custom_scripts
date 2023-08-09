#!/bin/bash
#SBATCH -J guppy
#SBATCH -A snic2022-22-85
#SBATCH -t 08:00:00
#SBATCH --exclusive
#SBATCH -p node
#SBATCH -N 1
#SBATCH -M snowy
#SBATCH --gpus=1
#SBATCH --gpus-per-node=1

# Script for basecalling on Nanopore long-read dRNA sequencing data using Guppy

module load bioinfo-tools
module load ont_fast5_api/3.1.6

guppy_basecaller="/domus/h1/amendez/ont-guppy-legacy/ont-guppy/bin/guppy_basecaller"

fast5_dir=$(realpath $1) # Peath to input raw Nanopore signal files in Fast5 format
outdir=$(realpath $2) # Output directory
flowcell="FLO-MIN106" # Sequencing flowcell
kit="SQK-RNA002" # Library preparation kit


# Running guppy basecaller using kit and flowcell information:
echo "Performing basecalling: "
${guppy_basecaller} -i ${fast5_dir} --save_path ${outdir} --flowcell ${flowcell} --kit ${kit} --recursive --reverse_sequence true --device "auto" --disable_pings --compress_fastq --verbose_logs

echo "Done."
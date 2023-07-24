#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 12
#SBATCH -t 20:00:00
#SBATCH -J xpdiff

# Script to perform m6A detection and differential modification rate analysis using xPore after Nanopore data preparation.
config_yaml=$(realpath $1) # Configuration file in YAML format (see https://xpore.readthedocs.io/en/latest/quickstart.html)
outdir=$(realpath $2) # Output directory

module load bioinfo-tools

echo "Activating xPore conda environment"

conda activate /proj/nb_project/private/conda_envs/nanopolish/

echo "Detecting m6A modification differences"

xpore diffmod --config ${config_yaml} \
--n_processes=${SLURM_NTASKS} \
--save_models=True

echo "Performing xPore's postprocessing command"

xpore postprocessing --diffmod_dir=${outdir}

echo "Done."
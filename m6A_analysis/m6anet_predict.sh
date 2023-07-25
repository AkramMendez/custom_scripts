#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 12
#SBATCH -t 08:00:00
#SBATCH -J m6aNetPred

dataprep_dir=$(realpath $1)
outdir=$(realpath $2)

module load bioinfo-tools

echo "Activating conda environment"

conda activate /proj/nb_project/private/conda_envs/nanopolish/

echo "Predicting m6A modification rates"

m6anet-run_inference --input_dir ${dataprep_dir} --out_dir ${outdir} --infer_mod-rate --n_processes ${SLURM_NTASKS}


echo "Done."
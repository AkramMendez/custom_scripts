#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 12
#SBATCH -t 08:00:00
#SBATCH -J m6aNetPrep

eventalign=$(realpath $1)
outdir=$(realpath $2)

module load bioinfo-tools

echo "Activating conda environment"

conda activate /proj/nb_project/private/conda_envs/nanopolish/


m6anet-dataprep --eventalign ${eventalign} \
                --out_dir ${outdir} \
                --n_processes ${SLURM_NTASKS}


echo "Done."

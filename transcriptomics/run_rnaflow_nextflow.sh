#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 16
#SBATCH -t 48:00:00
#SBATCH -J rnaflow

module load bioinfo-tools
#module load Nextflow
module load FastQC
module load fastp
#module load R_packages
module load samtools
module load subread
module load MultiQC

reads_file=$(realpath $1)
fastas_file=$(realpath $2)
gtf_file=$(realpath $3)
comparisons=$(realpath $4)
results_name=$5
outdir=$(realpath $6)

cd /crex/proj/nb_storage/private/rnaflow

#Note that NXF_HOME is set to $HOME/.nextflow
#Please change NXF_HOME to a place in your project directory (export NXF_HOME=yourprojectfolder)
export NXF_HOME=/crex/proj/nb_storage/private/rnaflow/nextflow_home

export CONDA_ENVS_PATH=/proj/nb_project/private/conda_envs
conda activate nextflow

echo "Running rnaflow"

nextflow run hoelzer-lab/rnaflow \
--reads ${reads_file} \
--genome ${fastas_file} \
--annotation ${gtf_file} \
-profile slurm,conda,latency \
--mode paired \
--strand 0 \
--skip_sortmerna \
--cores ${SLURM_NTASKS} \
--output ${results_name} \
--deg ${comparisons} \
-w ${outdir}

echo "Done."

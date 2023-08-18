#!/bin/bash
#SBATCH -J guppy
#SBATCH -A snic2020-15-304
#SBATCH -t 24:00:00
#SBATCH --exclusive
#SBATCH -p node
#SBATCH -N 1
#SBATCH -M snowy
#SBATCH --gres=gpu:1
#SBATCH --gpus-per-node=1
##for jobs shorter than 15 min (max 4 nodes):
#SBATCH --qos=short

#The time limit for jobs using GPUs is currently 3 days.

module load bioinfo-tools
module load ont_fast5_api/3.1.6

guppy_basecaller="/domus/h1/amendez/ont-guppy-legacy/ont-guppy/bin/guppy_basecaller"

fast5_dir=$(realpath $1)
outdir=$(realpath $2)
filename_base=$3
flowcel="FLO-MIN106"
kit="SQK-RNA002"

# Running basecaller using preloaded configuration files:
# ${guppy_basecaller} -i ${fast5_dir} --save_path ${outdir} -c ${config_file} --recursive --device "auto" --cpu_threads_per_caller 1 --gpu_runners_per_device 8 --num_callers ${SLURM_NTASKS} --disable_pings --compress_fastq --verbose_logs 

# Running guppy basecaller if using kit and flowcell information:
${guppy_basecaller} -i ${fast5_dir} --save_path ${outdir} --flowcell ${flowcell} --kit ${kit} --recursive --device "auto" --cpu_threads_per_caller 1 --gpu_runners_per_device 1 --num_callers ${SLURM_NTASKS} --disable_pings --compress_fastq --verbose_logs

#Without GPU usage:
${guppy_basecaller} -i ${fast5_dir} --save_path ${outdir} --flowcell ${flowcell} --kit ${kit} --recursive --disable_pings --compress_fastq --cpu_threads_per_caller 1 --num_callers ${SLURM_NTASKS} --verbose_logs


echo "Done."

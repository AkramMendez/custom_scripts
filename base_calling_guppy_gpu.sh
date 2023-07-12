#!/bin/bash
#SBATCH -J guppy
#SBATCH -A snic2020-15-304
#SBATCH -t 08:00:00
#SBATCH --exclusive
#SBATCH -p node
#SBATCH -N 1
#SBATCH -M snowy
#SBATCH --gpus=1
#SBATCH --gpus-per-node=1

#The time limit for jobs using GPUs is currently 3 days.

module load bioinfo-tools
module load ont_fast5_api/3.1.6

guppy_basecaller="/domus/h1/amendez/ont-guppy-legacy/ont-guppy/bin/guppy_basecaller"

fast5_dir=$(realpath $1)
outdir=$(realpath $2)
filename_base=""
#config_file=""
#multi_reads_dir=""
flowcell="FLO-MIN106"
kit="SQK-RNA002"

# Script for pre-processing Oxford Nanopore dRNA sequencing data using the Guppy basecaller.
# The Fast5 raw files can be processed either using GPUs or CPUs and using flowcell and kit model information to generate bascalled FASTQ files.

#Running basecaller using preloaded configuration files:

#${guppy_basecaller} -i ${fast5_dir} --save_path ${outdir} -c ${config_file} --recursive --device "auto" --cpu_threads_per_caller 1 --gpu_runners_per_device 1 --num_callers ${SLURM_NTASKS} --disable_pings --com$

# Running guppy basecaller if using kit and flowcell information:
${guppy_basecaller} -i ${fast5_dir} --save_path ${outdir} --flowcell ${flowcell} --kit ${kit} --recursive --reverse_sequence true --device "auto" --disable_pings --compress_fastq --verbose_logs

#${guppy_basecaller} -i ${fast5_dir} --save_path ${outdir} --flowcell ${flowcell} --kit ${kit} --recursive --disable_pings --compress_fastq --verbose_logs

#Note from:
#https://gist.github.com/sirselim/2ebe2807112fae93809aa18f096dbb94
# "* When performing GPU basecalling there is always one CPU support
#  thread per GPU caller, so the number of callers
#  (--num_callers) dictates the maximum number of CPU threads used.
#  * Max chunks per runner (--chunks_per_runner): The maximum number of
#  chunks which can be submitted to a single neural network runner
#  before it starts computation. Increasing this figure will increase
#  GPU basecalling performance when it is enabled.
#  * Number of GPU runners per device (--gpu_runners_per_device): The number of
#  neural network runners to create per CUDA device. Increasing this
#  number may improve performance on GPUs with a large number of
#  compute cores, but will increase GPU memory use. This option only

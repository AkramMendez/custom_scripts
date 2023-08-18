#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 16
#SBATCH -t 20:00:00
#SBATCH -J ribodetec

module load bioinfo-tools
module load TrimGalore/0.6.1

export CONDA_ENVS_PATH=/proj/nb_project/private/conda_envs
conda create -n ribodetector python=3.8
conda activate /proj/nb_project/private/conda_envs/ribodetector


inputdir=$(realpath $1)
outdir=$(realpath $2)


cd ${SNIC_TMP}

echo "Copying raw files to ${SNIC_TMP}"

rsync -avP ${inputdir}/ .

for i in $(ls | grep -oP "(.*)(?=_R(1|2)_001.fastq.gz)" | sort | uniq)
do 
sample=${i}; 
echo "Processing ${sample}:"
echo "Removing ribosomal reads..."

ribodetector_cpu -t ${SLURM_NTASKS} \
-l 100 \
-i ${sample}_R1_001.fastq.gz ${sample}/_R2_001.fastq.gz \
-e rrna \
--chunk_size 256 \
-o ${sample}_nonrrna_R1.fastq.gz ${sample}_nonrrna_R2.fastq.gz

echo "Trimming fastq files"

R1=$(echo "${sample}_nonrrna_R1.fastq.gz")
R2=$(echo "${sample}_nonrrna_R2.fastq.gz")

echo "Trimming files ${R1} and ${R2}"

trim_galore ${R1} ${R2} --paired --length 20 --cores ${SLURM_NTASKS} -o ${outdir}
done 

echo "Done."

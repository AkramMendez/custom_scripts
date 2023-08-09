#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 16
#SBATCH -t 12:00:00
#SBATCH -J hisat2_idx


module load bioinfo-tools
module load HISAT2/2.2.1

ref=$(realpath $1)
name=$2
outdir=$(realpath $3)

cd ${SNIC_TMP}
cp ${ref} ${SNIC_TMP}

genome=$(ls *.fa)
echo "Building Hisat2 index"


hisat2-build $genome $name -p ${SLURM_NTASKS}

echo "Listing files"
ls ${SNIC_TMP}

echo "Moving index files to ${outdir}"
mv *.ht2 ${outdir}

echo "Done"

#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 04:00:00
#SBATCH -J findMotifsHomer


module load bioinfo-tools
module load HOMER/4.11

ref=$(realpath $1) #reference fasta file
inputdir=$(realpath $2)
outdir=$(realpath $3)
sufif=$4 #Sufix to search directories within the main peak calling directory, each folder containing the peak calling results for a single sample
size=$5
len=$6



for i in $(find ${inputdir} -maxdepth 1 -type d -iname "${sufix}*")
do
name=$(basename ${i})

if [[ ! -d "${outdir}/${name}" ]]
then
    echo "Making output directory: ${outdir}/${name}"
    mkdir -p ${outdir}/${name}
fi

echo "Finding motifs for ${name}"

findMotifsGenome.pl ${i}/*.narrowPeak ${ref} ${outdir}/${name} \
-size ${size} \
-len ${len} \
-basic \
-rna \
-p ${SLURM_NTASKS}

done

echo "Done."
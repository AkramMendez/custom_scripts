#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 6
#SBATCH -t 06:00:00
#SBATCH -J hisat2_aln


#module load bioinfo-tools
#module load samtools
#module load HISAT2/2.2.1

R1=$1
sample_name=$2
outdir=$3


n=$(echo "$(($SLURM_NTASKS-1))")

hisat2_index_basename=/proj/nb_project/private/genomes/hisat2_indexes/mm10_hisat2/mm10/genome

#echo "Creating output directory: ${outdir}"
#mkdir ${outdir}

echo "Mapping reads to index genome using HISAT2"

hisat2 -p ${n} -x ${hisat2_index_basename} -U ${R1} --rna-strandness F --summary-file ${outdir}/hisat2_summary_${sample_name}.txt | samtools view -@ ${n} -bS - > ${outdir}/${sample_name}.hisat2.bam

echo "Mapping finished."

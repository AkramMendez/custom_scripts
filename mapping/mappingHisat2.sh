#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 06:00:00
#SBATCH -J hisat2_aln


module load bioinfo-tools
module load HISAT2/2.2.1

R1=$1
R2=$2
sample_name=$3
outdir=$4


n=$(echo "$(($SLURM_NTASKS-1))")


hisat2_index_basename=/proj/nb_project/private/hisat2indexes/hisat2index_chlSab2_wuhCor1_ecoliK12/chlSab2_wuhCor1_ecoliK12

#echo "Creating output directory: ${outdir}"
#mkdir ${outdir}

echo "Mapping reads to index genome using HISAT2"

hisat2 -p ${n} -x ${hisat2_index_basename} -1 ${R1} -2 ${R2} --rna-strandness RF --summary-file ${outdir}/hisat2_summary_${sample_name}.txt | samtools view -h -bS - > ${outdir}/${sample_name}.hisat2.bam

samtools sort ${outdir}/${sample_name}.hisat2.bam -o ${outdir}/${sample_name}.hisat2.sorted.bam -@ ${SLURM_NTASKS}

samtools index ${outdir}/${sample_name}.hisat2.sorted.bam -@ ${SLURM_NTASKS}

echo "Mapping finished. BAM mapping results saved in ${outdir}"

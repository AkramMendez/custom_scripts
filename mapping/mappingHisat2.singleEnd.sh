#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 6
#SBATCH -t 06:00:00
#SBATCH -J hisat2_aln


module load bioinfo-tools
module load HISAT2/2.2.1
module load samtools


R1=$(realpath $1)
sample_name=$(basename $R1)
hisat2_index_basename=$(realpath $3)
outdir=$(realpath $4)


n=$(echo "$(($SLURM_NTASKS-1))")



#echo "Creating output directory: ${outdir}"
#mkdir ${outdir}

cd ${SNIC_TMP}

echo "Mapping reads to index genome using HISAT2"

hisat2 -p ${n} -x ${hisat2_index_basename} -U ${R1} --rna-strandness R --summary-file hisat2_summary_${sample_name}.txt | samtools view -@ ${n} -bS - | samtools sort -o ${sample_name}.hisat2.sorted.bam -@ ${n}

echo "Mapping finished."

echo "Listing files"

ls ${SNIC_TMP}

echo "Copying files to ${outdir}"

cp ${SNIC_TMP}/*.bam ${outdir}
cp ${SNIC_TMP}/*.txt ${outdir}

echo "Done."
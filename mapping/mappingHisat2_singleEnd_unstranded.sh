#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 6
#SBATCH -t 06:00:00
#SBATCH -J hisat2_aln


module load bioinfo-tools
module load HISAT2/2.2.1
module load samtools


inputdir=$(realpath $1)
hisat2_index_basename=$(realpath $2)
outdir=$(realpath $3)


n=$(echo "$(($SLURM_NTASKS-1))")


cd ${SNIC_TMP}

cp ${inputdir}/*.fastq.gz .

for i in $(ls *.fastq.gz)
do
	R1=${i}
	sample_name=$(${i}%.fastq.gz)
	echo "Mapping reads to index genome using HISAT2"

	hisat2 -p ${n} -x ${hisat2_index_basename} -U ${R1} --summary-file hisat2_summary_${sample_name}.txt | samtools view -@ ${n} -bS - | samtools sort -o ${sample_name}.hisat2.sorted.bam -@ ${n}

	samtools index ${${sample_name}.hisat2.sorted.bam} -@ ${SLURM_NTASKS}	
done

echo "Mapping finished."
echo "Listing files"

ls ${SNIC_TMP}

echo "Copying files to ${outdir}"

cp ${SNIC_TMP}/*.bam ${outdir}
cp ${SNIC_TMP}/*.txt ${outdir}

echo "Done."
#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
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

cp ${inputdir}/fq.gz .

for R1 in $(ls *.fq.gz)
do
	sample_name=$(basename $R1)

	echo "Mapping reads to index genome using HISAT2"

	hisat2 -p ${n} -x ${hisat2_index_basename} -U ${R1} --rna-strandness R --summary-file hisat2_summary_${sample_name}.txt | samtools view -@ ${n} -bS - | samtools sort -o ${sample_name}.hisat2.sorted.bam -@ ${n}
	
	echo "Mapping finished."
	echo "Listing files"

	ls ${SNIC_TMP}

	echo "Copying files to ${outdir}"
 
	cp ${SNIC_TMP}/*.bam ${outdir}
	cp ${SNIC_TMP}/*.txt ${outdir}

done


echo "Done."

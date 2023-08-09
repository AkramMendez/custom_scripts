#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 18:00:00
#SBATCH -J hisat2_aln_markdup

module load bioinfo-tools
module load HISAT2/2.2.1
module load samtools
module load picard/2.23.4


inputdir=$(realpath $1)
hisat2_index_basename=$(realpath $2)
outdir=$(realpath $3)


n=$(echo "$(($SLURM_NTASKS-1))")



#echo "Creating output directory: ${outdir}"
#mkdir ${outdir}

cd ${SNIC_TMP}
cp ${inputdir}/*.fq.gz .

for R1 in $(ls *.fq.gz)
do
	sample_name=$(basename $R1)
	echo "Mapping reads to index genome using HISAT2"
	hisat2 -p ${n} -x ${hisat2_index_basename} -U ${R1} --rna-strandness F --summary-file hisat2_summary_${sample_name%.fq.gz}.txt | samtools view -@ ${n} -bS - | samtools sort -o ${sample_name%.fq.gz}.hisat2.sorted.bam -@ ${n}
	
done

echo "Mapping finished."
	echo "Listing files"

for bam in $(ls *.sorted.bam)
do 
	echo "Marking duplicates, file: ${bam}"

	java -jar $PICARD_ROOT/picard.jar MarkDuplicates I=${bam} O=${bam%.sorted.bam}.markdup.bam M=${bam%.sorted.bam}.markdup.metrics.txt TAGGING_POLICY=All

	echo "Saving marked alignment to ${outdir}"

	cp ${bam%.sorted.bam}.markdup.bam ${outdir}
	cp ${bam%.sorted.bam}.markdup.metrics.txt ${outdir}

done


	ls ${SNIC_TMP}
	echo "Copying files to ${outdir}"

	cp ${SNIC_TMP}/*.markdup.bam ${outdir}
	cp ${SNIC_TMP}/*.txt ${outdir}

echo "Done."
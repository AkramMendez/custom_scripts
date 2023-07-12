#!/bin/bash -l
#SBATCH -A naiss2023-22-153
#SBATCH -p core -n 8
#SBATCH -t 24:00:00
#SBATCH -J markDuplicates

module load bioinfo-tools
module load picard/2.23.4
module load samtools/1.12

inputdir=$(realpath $1)  # Path to original sorted BAM files (can be coordinate-sorted or query-sorted)
outdir=$(realpath $2) # Output directory

# Script for marking duplicated reads in BAM files using Picard (https://gatk.broadinstitute.org/hc/en-us/articles/360037052812-MarkDuplicates-Picard-)

cd $SNIC_TMP

cp ${inputdir}/*.sorted.bam .

for bam in $(ls *.sorted.bam)
do 
	echo "Marking duplicates, file: ${bam}"

	java -jar $PICARD_ROOT/picard.jar MarkDuplicates I=${bam} O=${bam%.sorted.bam}.markdup.bam M=${bam%.sorted.bam}.markdup.metrics.txt TAGGING_POLICY=All

	echo "Saving marked alignment to ${outdir}"

	cp ${bam%.sorted.bam}.markdup.bam ${outdir}
	cp ${bam%.sorted.bam}.markdup.metrics.txt ${outdir}

done

echo "Done."
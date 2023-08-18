#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 24:00:00
#SBATCH -J markDuplicates

# Script for marking duplicates from BAM files using Picard.
module load bioinfo-tools
module load picard/2.23.4
module load samtools/1.12

inputdir=$(realpath $1) # Input directory with sorted BAM files
outdir=$(realpath $2) # Output direcory

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
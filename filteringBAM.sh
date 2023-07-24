#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 03:00:00
#SBATCH -J filtering_bam

# Script for filtering BAM files to remove duplicate and unmapped reads.
module load bioinfo-tools
module load sambamba/0.7.1
module load samtools/1.12

inputdir=$(realpath $1) # Input directory path with duplicate-marked BAM files
outdir=$(realpath $2) # Output directory

cd $SNIC_TMP

cp ${inputdir}/*.markdup.bam .

for bam in $(ls *.markdup.bam)
do 
	echo "Filtering aln file: ${bam}"

	sambamba view -q -h -f bam -t ${SLURM_NTASKS} -F "not unmapped and not duplicate" ${bam} > ${bam%.markdup.bam}.mapped.nodup.bam

	echo "Sorting and indexing filtered alignment"
	samtools sort ${bam%.markdup.bam}.mapped.nodup.bam -o ${bam%.markdup.bam}.mapped.nodup.sorted.bam -@ ${SLURM_NTASKS}

	samtools index ${bam%.markdup.bam}.mapped.nodup.sorted.bam -@ ${SLURM_NTASKS}

	echo "Saving filtered alignment to ${outdir}"
	cp ${bam%.markdup.bam}.mapped.nodup.sorted.bam ${outdir}
	cp *.bai ${outdir}

done

echo "Done."
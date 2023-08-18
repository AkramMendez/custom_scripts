#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 03:00:00
#SBATCH -J split_bam

# Script to split BAM files by strand
module load bioinfo-tools
module load sambamba/0.7.1
module load samtools/1.12

inputdir=$(realpath $1) # Input directory with unsplit BAM files
outdir=$(realpath $2) # Output directory

cd $SNIC_TMP

cp ${inputdir}/*.bam .

for bam in $(ls *.bam)
do 
	echo "Filtering aln file: ${bam}"

	sambamba view -q -h -f bam -t ${SLURM_NTASKS} -F "[XS]=='+'" ${bam} | samtools sort -o ${bam%.bam}.fwd.bam -@ ${SLURM_NTASKS}
	sambamba view -q -h -f bam -t ${SLURM_NTASKS} -F "[XS]=='-'" ${bam} | samtools sort -o ${bam%.bam}.rev.bam -@ ${SLURM_NTASKS}

	echo "Indexing filtered alignment"
	samtools index ${bam%.bam}.fwd.bam -@ ${SLURM_NTASKS}
	samtools index ${bam%.bam}.rev.bam -@ ${SLURM_NTASKS}

	echo "Saving filtered alignment to ${outdir}"
	cp *.fwd.bam ${outdir}
	cp *.rev.bam ${outdir}
	cp *.bai ${outdir}
	
done

echo "Done."
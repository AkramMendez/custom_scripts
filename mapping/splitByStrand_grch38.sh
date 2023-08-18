#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 03:00:00
#SBATCH -J split_bam

module load bioinfo-tools
module load sambamba/0.7.1
module load samtools/1.12

inputdir=$(realpath $1)
outdir=$(realpath $2)

cd $SNIC_TMP

cp ${inputdir}/*.bam .

for bam in $(ls *.bam)
do 
	echo "Filtering aln file: ${bam}"

	sambamba view -q -h -f bam -t ${SLURM_NTASKS} -F "[XS]=='+' and not unmapped and not duplicate and [NH]==1" ${bam} | samtools sort -o ${bam%.bam}.fwd.bam -@ ${SLURM_NTASKS}
	sambamba view -q -h -f bam -t ${SLURM_NTASKS} -F "[XS]=='-' and not unmapped and not duplicate and [NH]==1" ${bam} | samtools sort -o ${bam%.bam}.rev.bam -@ ${SLURM_NTASKS}

	echo "Indexing filtered alignment"
	samtools index ${bam%.bam}.fwd.bam -@ ${SLURM_NTASKS}
	samtools index ${bam%.bam}.rev.bam -@ ${SLURM_NTASKS}

	echo "Saving filtered alignment to ${outdir}"
	cp *.fwd.bam ${outdir}
	cp *.rev.bam ${outdir}
	cp *.bai ${outdir}
	
done

echo "Done."
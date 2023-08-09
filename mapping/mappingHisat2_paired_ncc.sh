#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 24:00:00
#SBATCH -J hisat2_aln

module load bioinfo-tools
module load HISAT2/2.2.1
module load samtools
module load picard/2.23.4
module load sambamba


inputdir=$(realpath $1)
hisat2_index_basename=$(realpath $2)
outdir=$(realpath $3)


n=$(echo "$(($SLURM_NTASKS-1))")


#echo "Creating output directory: ${outdir}"
#mkdir ${outdir}

cd ${SNIC_TMP}
cp ${inputdir}/*.fq.gz .

echo "Mapping reads"

for read in $(ls | grep -oP "(.*)(?=_\d+.fq.gz)" | sort | uniq)
do
	sample_name=$(basename $read)
	R1=$(echo "${read}"_1.fq.gz)
	R2=$(echo "${read}"_2.fq.gz)
	
	echo "Mapping reads to index genome using HISAT2"
	# Parameters reported in BGI's page: --sensitive --no-discordant --no-mixed -I 1 -X 1000
	#hisat2 -p ${n} -x ${hisat2_index_basename} -U ${R1} --summary-file hisat2_summary_${sample_name%.fq.gz}.txt | samtools view -@ ${n} -bS - | samtools sort -o ${sample_name%.fq.gz}.hisat2.sorted.bam -@ ${n}
	hisat2 -p ${n} -x ${hisat2_index_basename} --sensitive --no-discordant --no-mixed -I 1 -X 1000 -1 ${R1} -2 ${R2} --summary-file hisat2_summary_${sample_name%.fq.gz}.txt | samtools view -@ ${n} -bS - | samtools sort -o ${sample_name%.fq.gz}.hisat2.sorted.bam -@ ${n}

	samtools index ${sample_name%.fq.gz}.hisat2.sorted.bam -@ ${SLURM_NTASKS}

	echo "Copying hisat2 summary: hisat2_summary_${sample_name%.fq.gz}.txt  to ${outdir}"
	
	cp hisat2_summary_${sample_name%.fq.gz}.txt ${outdir}

done

echo "Mapping finished"

echo "Marking duplicates"

for bam in $(ls *.sorted.bam)
do 
	echo "Marking duplicates, file: ${bam}"

	java -jar $PICARD_ROOT/picard.jar MarkDuplicates I=${bam} O=${bam%.sorted.bam}.markdup.bam M=${bam%.sorted.bam}.markdup.metrics.txt TAGGING_POLICY=All

	
	echo "Indexing marked-duplicates alignment:"

	samtools index ${bam%.sorted.bam}.markdup.bam -@ ${SLURM_NTASKS}

	#echo "Saving marked alignment to ${outdir}"

	#cp ${bam%.sorted.bam}.markdup.bam ${outdir}

	#cp ${bam%.sorted.bam}.markdup.bam.bai ${outdir}

	echo "Copying ${bam%.sorted.bam}.markdup.metrics.txt to ${outdir}"

	cp ${bam%.sorted.bam}.markdup.metrics.txt ${outdir}

done

echo "Filtering alignments"

for markbam in $(ls *.markdup.bam)
do 
	echo "Filtering aln file: ${markbam}"

	sambamba view -q -h -f bam -t ${SLURM_NTASKS} -F "[NH]==1 and not unmapped and not duplicate" ${markbam} > ${markbam%.markdup.bam}.nodup.uniq.bam

	echo "Sorting and indexing filtered alignment"
	samtools sort ${markbam%.markdup.bam}.nodup.uniq.bam -o ${markbam%.markdup.bam}.nodup.uniq.sorted.bam -@ ${SLURM_NTASKS}

	samtools index ${markbam%.markdup.bam}.nodup.uniq.sorted.bam -@ ${SLURM_NTASKS}

	echo "Saving filtered alignment to ${outdir}"
	cp ${markbam%.markdup.bam}.nodup.uniq.sorted.bam ${outdir}
	cp ${markbam%.markdup.bam}.nodup.uniq.sorted.bam.bai ${outdir}

done


	#ls ${SNIC_TMP}
	#echo "Copying files to ${outdir}"

	#cp ${SNIC_TMP}/*.markdup.bam ${outdir}
	#cp ${SNIC_TMP}/hisat2_summary_*.txt ${outdir}

echo "Done."
#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 06:00:00
#SBATCH -J markdup

module load bioinfo-tools
module load samtools
module load sambamba
module load bamtools

#sam=$1
bam=$1
#n=$2
#outdir=$(realpath $2)

#bam=$(basename ${sam/.sam/.bam})
n=${SLURM_NTASKS}
echo "Converting SAM to BAM"


#samtools view -bS ${sam} -@ ${n} -o ${outdir}/${bam}

#cd ${outdir}

echo "Collate bam by name"
samtools collate -o ${bam%.bam}.collate.bam ${bam} --threads ${n}

echo "Running fixmante"
samtools fixmate -m ${bam%.bam}.collate.bam ${bam%.bam}.collate.fixmate.bam --threads ${n}

echo "Sorting fixmate bam file"
samtools sort -o ${bam%.bam}.collate.fixmate.sorted.bam ${bam%.bam}.collate.fixmate.bam --threads ${n}

echo "Marking duplicates"
samtools markdup ${bam%.bam}.collate.fixmate.sorted.bam ${bam%.bam}.collate.fixmate.sorted.markdup.bam --threads ${n}

echo "Getting stats from markdup bam file"

bamtools stats -in ${bam%.bam}.collate.fixmate.sorted.markdup.bam > bamtools.stats.${bam%.bam}.markdup.bam.txt

echo "Removing duplicates"

sambamba view -h -t ${n} -f bam -F "not unmapped and not duplicate" ${bam%.bam}.collate.fixmate.sorted.markdup.bam > ${bam%.bam}.nodup.bam

#echo "Split by strand"
# Note: When using Hisat2, the mapping file includes NH and XS flags, the first indicating the number of hits to one or multiple locations, 
# and the XS indicating the strant to where the read mapped to ('+' or '-'). You can filter out multimapping reads by selecting NH ==1, or allow number of multimapped reads NH>1.
#sambamba view -h -t ${n} -f bam -F "[NH]==1 and not unmapped and not duplicate" ${bam%.bam}.collate.fixmate.sorted.markdup.bam > ${bam%.bam}.multihit1.nodup.bam
#sambamba view -h -t ${n} -f bam -F "[XS]== '-' and not unmapped and not duplicate" ${bam%.bam}.collate.fixmate.sorted.markdup.bam > ${bam%.bam}.nodup.rev.bam
#sambamba view -h -t ${n} -f bam -F "[XS]== '+' and not unmapped and not duplicate" ${bam%.bam}.collate.fixmate.sorted.markdup.bam > ${bam%.bam}.nodup.fwd.bam

echo "Split by strand"

# Filtering by strand (for Paired-end data):
#samtools view -h ${bam%.bam}.nodup.bam -@ ${n} | awk '{if($1~/^@/ || $2 ~ /^(16|83|163)$/){print}}' | samtools view -h -bS -@ ${n} - > ${bam%.bam}.nodup.posfwd.bam

#samtools view -h ${bam%.bam}.nodup.bam -@ ${n} | awk '{if($1~/^@/ || $2 ~ /^(0|99|147)$/){print}}' | samtools view -h -bS -@ ${n} - > ${bam%.bam}.nodup.negrev.bam

echo "Sorting BAM files"
samtools sort ${bam%.bam}.nodup.posfwd.bam -o ${bam%.bam}.nodup.sorted.posfwd.bam --threads ${n}

samtools sort ${bam%.bam}.nodup.negrev.bam -o ${bam%.bam}.nodup.sorted.negrev.bam --threads ${n}

echo "Indexing BAM files"

samtools index ${bam%.bam}.nodup.sorted.posfwd.bam -@ ${n}

samtools index ${bam%.bam}.nodup.sorted.negrev.bam -@ ${n}

echo "Removing intermediate BAM files"
#rm ${sam}
#rm ${bam}
rm ${bam%.bam}.collate.bam
rm ${bam%.bam}.collate.fixmate.bam
rm ${bam%.bam}.collate.fixmate.sorted.bam
#rm ${bam%.bam}.collate.fixmate.sorted.markdup.bam
rm ${bam%.bam}.nodup.bam
rm ${bam%.bam}.nodup.posfwd.bam
rm ${bam%.bam}.nodup.negrev.bam

echo "Done."

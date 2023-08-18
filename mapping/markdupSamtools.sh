#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 06:00:00
#SBATCH -J markdup

module load bioinfo-tools
module load samtools
module load sambamba

sam=$1
outdir=$(realpath $2)

bam=${sam%.sam}.bam

echo "Converting SAM to BAM"


samtools view -S -b ${sam} --threads ${n} > ${bam}

echo "Collate bam by name"
samtools collate -o ${bam%.bam}.collate.bam ${bam} --threads ${n}

echo "Running fixmante"
samtools fixmate -m ${bam%.bam}.collate.bam ${bam%.bam}.collate.fixmate.bam --threads ${n}

echo "Sorting fixmate bam file"
samtools sort -o ${bam%.bam}.collate.fixmate.sorted.bam ${bam%.bam}.collate.fixmate.bam --threads ${n}

echo "Marking duplicates"
samtools markdup ${bam%.bam}.collate.fixmate.sorted.bam ${bam%.bam}.collate.fixmate.sorted.markdup.bam --threads ${n}

echo "Removing duplicates"

sambamba view -h -t ${n} -f bam -F "not null and not unmapped and not duplicate" ${bam%.bam}.collate.fixmate.sorted.markdup.bam > ${bam%.bam}.nodup.bam
#sambamba view -h -t ${n} -f bam -F "[NH]==1 and not unmapped and not duplicate" ${bam%.bam}.collate.fixmate.sorted.markdup.bam > ${bam%.bam}.multihit1.nodup.bam

samtools sort -o ${bam%.bam}.nodup.sorted.bam ${bam%.bam}.nodup.bam --threads ${n}

echo "Removing intermediate BAM files"
rm ${sam}
rm ${bam}
rm ${bam%.bam}.collate.bam
rm ${bam%.bam}.collate.fixmate.bam
rm ${bam%.bam}.collate.fixmate.sorted.bam
rm ${bam%.bam}.collate.fixmate.sorted.markdup.bam
rm ${bam%.bam}.nodup.bam
echo "Done."

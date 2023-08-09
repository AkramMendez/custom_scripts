#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 06:00:00
#SBATCH -J star_mapping

module load bioinfo-tools
module load gcc/10.2.0 star/2.7.2b

IDX=${1}
READS=${2}
OUT=${3}


for READ in $(ls ${READS/*.fq.gz})
do
	PREF=echo "${OUT}/${read%.fq.gz}_"
	STAR --genomeDir ${IDX} \
	--runThreadN 10 \
	--readFilesIn ${READ} \
	--readFilesCommand zcat \
	--outFileNamePrefix ${PREF} \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes Standard
done

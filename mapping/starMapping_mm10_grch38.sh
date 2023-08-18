#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 10
#SBATCH -t 24:00:00
#SBATCH -J STAR

module load bioinfo-tools
module load gcc star/2.7.9a

idx=$(realpath $1)
inputdir=$(realpath $2)
outdir=$(realpath $3)
species=$4

#Preloaded indexes in Uppmax:
#/sw/data/igenomes/Mus_musculus/UCSC/mm10/Sequence/STARIndex
#/sw/data/igenomes/Homo_sapiens/UCSC/hg38/Sequence/STARIndex

cd ${SNIC_TMP}
echo "Copying files to ${SNIC_TMP}"

cp ${inputdir}/*fq.gz .


samples=$(ls *.fq.gz | sed -E 's/_(1|2).fq.gz//g' | sort | uniq | tr "\n" " ")


for sample in ${samples[@]}
do
	prefix=$(echo "${sample}_${species}")
	STAR --genomeDir ${idx} \
	--runThreadN ${SLURM_NTASKS} \
	--readFilesIn ${sample}_1.fq.gz ${sample}_2.fq.gz \
	--readFilesCommand zcat \
	--outFileNamePrefix ${prefix} \
	--outSAMtype BAM SortedByCoordinate \
	--outSAMunmapped Within \
	--outSAMattributes NM #Set NM attribute to get extra information if you plan to use XenofilteR: https://github.com/NKI-GCF/XenofilteR
done


echo "Copying BAM files to ${outdir}"
cp *.bam ${outdir}
cp *.bai ${outdir}

echo "Done."
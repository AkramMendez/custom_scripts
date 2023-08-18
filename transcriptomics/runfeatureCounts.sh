#!/bin/bash -l
#SBATCH -A snic2022-2-85
#SBATCH -p core -n 8
#SBATCH -t 08:00:00
#SBATCH -J featureCounts

module load bioinfo-tools
module load subread/2.0.0


gtf=$(realpath $1)
inputdir=$(realpath $2)
outdir=$(realpath $3)
name=$4

cd $SNIC_TMP

cp ${inputdir}/*.bam .
cp ${inputdir}/*.bai .


date_stamp=$(date +"%y-%m-%d")

featureCounts -t exon -g gene_id \
-s 1 \
-a ${gtf} \
-o featureCounts.${name}_${date_stamp}.txt \
$(ls ${inputdir}/*.bam | tr "\n" " ") \
-T ${SLURM_NTASKS} \
2> ${name}.featureCounts.${name}_stranded_${date_stamp}.log


echo "Moving files to ${outdir}"

cp *.txt ${outdir}
cp *.log ${outdir}

echo "Done."
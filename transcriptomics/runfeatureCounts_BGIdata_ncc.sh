#!/bin/bash -l
#SBATCH -A snic2020-15-304
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

#-s 0, for BGI dataset which is unstranded
#s 1, for stranded data, if m6A SMARTer sequences were reverse complemented use this, otherwise use s 2 (reverse stranded) for the SMARTer Takara Pico Inut v2 reads aligned in the original orientation.
featureCounts -t exon -g gene_id \
-s 0 \
-a ${gtf} \
-p \
-B \
-o featureCounts.${name}_${date_stamp}.txt \
$(ls *.bam | tr "\n" " ") \
-T ${SLURM_NTASKS} \
--verbose \
2> ${name}.featureCounts.${name}_stranded_${date_stamp}.log


echo "Moving files to ${outdir}"

cp *.txt ${outdir}
cp *.log ${outdir}

echo "Done."


#featureCounts -t exon -g gene_id -s 1 -a ${gtf} -p -B -o featureCounts.${name}_${date_stamp}.txt $(ls *.bam | tr "\n" " ") -T ${SLURM_NTASKS} --verbose 2> ${name}.featureCounts.${name}_stranded_${date_stamp}.log

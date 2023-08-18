#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 28:00:00
#SBATCH -J revcomp_seqkit

module load bioinfo-tools
module load SeqKit/0.15.0


inputdir=$(realpath $1)
outdir=$(realpath $2)

cd $SNIC_TMP
cp ${inputdir}/*.gz .

for sample in $(ls *.gz)
do
	echo "Processing sample ${sample}"

	seqkit seq --seq-type DNA -j ${SLURM_NTASKS} -r -p ${sample} | gzip > ${sample%.fq.gz}.revcomp.fq.gz
	echo "Copying revcomp file to ${outdir}:"
	cp ${sample%.fq.gz}.revcomp.fq.gz ${outdir}
done


echo "Done." 
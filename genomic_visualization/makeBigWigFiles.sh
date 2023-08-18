#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 08:00:00
#SBATCH -J bamCoverage

module load bioinfo-tools
module load deepTools

inputdir=$(realpath $1)
outdir=$(realpath $2)

cd $SNIC_TMP

cp ${inputdir}/*.bam .
cp ${inputdir}/*.bai .

bamfiles=$(ls *.bam)

n=${SLURM_NTASKS}

for bam in ${bamfiles}
do
	echo "Making BigWig for ${bam}"
	outbam=${bam%.bam}.bs10.bw
	bamCoverage -b ${bam} -p ${n} --binSize 10 --outFileFormat bigwig -o ${outbam}
done

echo "Copying files to ${outdir}"

cp *.bw ${outdir}

echo "Done."
#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 6
#SBATCH -t 06:00:00
#SBATCH -J bamCoverage

module load bioinfo-tools
module load deepTools

samplesdir=$1
outdir=$2

n=$(echo "$(($SLURM_NTASKS-1))")

plotBigWig(){

echo "Generating coverage bedgraph for sample: ${sample}"
bamCoverage --bam ${bam} \
--outFileName ${outdir}/${bam%.bam}.bw \
--outFileFormat bigwig \
--normalizeUsing RPKM \
--binSize 1 \
--numberOfProcessors ${n} \
--verbose

}


for bam in $(ls ${samplesdir}/*.bam)
do
	plotBigWig ${bam} ${outdir} ${n}
done

echo "Done."
#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 16
#SBATCH -t 08:00:00
#SBATCH -J plotCor

module load bioinfo-tools
module load picard/2.23.4
module load sambamba/0.7.1
module load samtools/1.12
module load subread/2.0.0
module load deepTools

inputdir=$(realpath $1)
outdir=$(realpath $2)

cd ${SNIC_TMP}

echo "Copying fils to ${SNIC_TMP}"
cp ${inputdir}/*.nodup.uniq.bam* .

multiBamSummary bins --bamfiles $(ls *.nodup.uniq.bam | tr "\n" " ") \
--labels $(ls *.nodup.uniq.bam | sed 's/.hg38.cov2.nodup.uniq.bam//g' | tr "\n" " ") \
--outFileName multiBamSummary.allSamples.nodup.uniq.npz \
-p ${SLURM_NTASKS}

plotCorrelation --corData multiBamSummary.allSamples.nodup.uniq.npz \
--corMethod pearson \
--skipZeros \
--whatToPlot heatmap \
--plotFileFormat pdf \
--plotFile correlationPlot.allSamples.nodup.uniq.pdf \
--colorMap Blues


echo "Copying files to ${outdir}"

cp *.pdf ${outdir}
cp *.npz ${outdir}

#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 06:00:00
#SBATCH -J downsampPicard

# Script to perform downsampling on BAM files using Picard and provided scaling factors.
module load bioinfo-tools
module load picard/2.23.4
module load samtools/1.12

bam_file=$(realpath $1) # Input BAM file path
outdir=$(realpath $2) # Output directory
scaling_factor=$3 # Scaling factor for downsampling

cd $SNIC_TMP

cp ${bam_file} .

downsampleBAM(){
	bam=$1
	scalefactor=$2
	outdir=$3
	echo "Downsampling file: ${bam}"

	java -jar $PICARD_ROOT/picard.jar DownsampleSam -I ${bam} \
	-O ${outdir}/${bam%.bam.sorted}.spikeInNorm.bam \
	-STRATEGY HighAccuracy \
	-P ${scalefactor}	
}

downsampleBAM ${bam_file} ${scaling_factor} ${outdir}

echo "Done."

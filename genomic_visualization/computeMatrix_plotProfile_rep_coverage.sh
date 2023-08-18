#!/bin/bash -l
#SBATCH -A snic2022-15-85
#SBATCH -p core -n 10
#SBATCH -t 24:00:00
#SBATCH -J computeMat


module load bioinfo-tools
module load deepTools/3.3.2

inputdir=$(realpath $1) # Path to input bigwig track files
outdir=$(realpath $2) # Output directory
ref=$(realpath $3) # Reference GTF annotation file
name=$4 # Custom label for naming the output files

# Script for computing the coverage across gene feature regions +- 3kb and plotting the scaled-regions profiles from coverage bigwig tracks

cd ${SNIC_TMP}

cp ${inputdir}/*.bw .

computeMatrix scale-regions -S $(ls *.bw | tr "\n" " ") -R ${ref} -a 3000 -b 3000 --regionBodyLength 3000 --outFileName scale_regions_${name}.mat.gz --outFileNameMatrix scale_regions_${name}.tabular --samplesLabel $(ls *.bw | sed -E 's/.bw//g' | tr "\n" " ") --verbose -p ${SLURM_NTASKS} --skipZeros

plotProfile --matrixFile scale_regions_${name}.mat.gz --perGroup --plotFileFormat pdf --outFileName profile_scaleRegions_${name}.pdf --plotTitle "" --regionsLabel "" --outFileNameData profile_scaleRegions_${name}.tsv

echo "Copying result files to ${outdir}:"

cp *.gz ${outdir}
cp *.tabular ${outdir}
cp *.tsv ${outdir}
cp *.pdf ${outdir}

echo "Done."


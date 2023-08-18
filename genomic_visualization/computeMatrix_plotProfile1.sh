#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 12
#SBATCH -t 12:00:00
#SBATCH -J bamCoverage


module load bioinfo-tools
module load deepTools/3.3.2

ref=$(realpath $1) # Path to GTF reference annotation file (i.e gencode.v36.chr_patch_hapl_scaff.basic.annotation.gtf)
inputdir=$(realpath $2) # Path to bigwig track files
outdir=$(realpath $3) # Output directory

# Script for a computing the coverage across annotated genomic features using deepTools computeMatrix, scaled-regions coverage across gene annotations +- 3kb and plotted in PDF format.

cd $SNIC_TMP

cp ${inputdir}/*.bw .

echo "ComputeMatrix calculation"
computeMatrix scale-regions -S $(ls *.bw | tr "\n" " ") -R ${ref} \
-a 3000 \
-b 3000 \
--regionBodyLength 3000 \
--outFileName scale_regions.mat.gz \
--outFileNameMatrix scale_regions.tabular \
--samplesLabel $(ls *.bw | sed -E 's/.bw//g' | tr "\n" " ") \
--verbose \
-p ${SLURM_NTASKS} \
--skipZeros

echo "Computing plotProfile"
plotProfile --matrixFile scale_regions.mat.gz \
--perGroup \
--plotFileFormat pdf \
--outFileName ${outdir}/profile_scaleRegions.pdf \
--plotTitle "" \
--regionsLabel "" \
--outFileNameData ${outdir}/scale_regions.tsv

echo "Copying result files to ${outdir}:"

cp *.gz ${outdir}
cp *.tabular ${outdir}

echo "Done."

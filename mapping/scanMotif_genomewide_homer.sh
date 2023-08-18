#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 06:00:00
#SBATCH -J scanMotif



module load bioinfo-tools
module load HOMER


motif=$1
mismatch=$2
genome=$(realpath $3)
outdir=$(realpath $4)

echo "Scanning for motif genome-wide:"
name=${motif}

seq2profile.pl ${motif} ${mismatch} ${name} > ${name}.motif 

motif_file=${name}.motif

scanMotifGenomeWide.pl ${motif_file} ${genome} -bed > ${name}_homer_motifScanGenome.bed

echo "Copying motif file and scan results to ${outdir}"

cp ${name}_homer_motifScanGenome.bed ${outdir}
cp ${name}.motif ${outdir}

echo "Done."
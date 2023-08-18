#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 08:00:00
#SBATCH -J splitBam

module load bioinfo-tools
module load picard/2.23.4
module load samtools/1.12

inputdir=$(realpath $1) # Path to input unsplit BAM files
outdir=$(realpath $2) # Output directory

# Script for splitting a BAM file into smaller chunks using Picard. Useful for processing huge BAM files for downstream applications with a limited number of reads to process.
# The split alignments are named with the suffix 001,002 depending on the number of chunks defined by the SPLIT_TO_N_FILES parameter

cd $SNIC_TMP

cp ${inputdir}/*.bam .
cp ${inputdir}/*.bai .

echo "Splitting BAM files"
for bam in $(ls *.bam)
do
    java -jar $PICARD_ROOT/picard.jar SplitSamByNumberOfReads \
    I=${bam} \
    OUTPUT=${outdir} \
    SPLIT_TO_N_FILES=2
done

echo "Done."

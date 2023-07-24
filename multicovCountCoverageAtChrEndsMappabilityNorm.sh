#!/usr/bin/bash
# Script to calculate the number of reads inside the chromosome end coordinates for each alignment provided in BED format using multicov from BEDtools.
# The counts are further used to calculate the ammount of expression from each telomere end following mappability and library size normalization.

input_files=$(realpath $1) # Input BAM files
chrEnds_P=$(realpath $2) # BED file with chromosome end (P arm) coordinates
chrEnds_Q=$(realpath $3) # BED file with chromosome end (Q arm) coordinates
outdir=$(realpath $4) # Output directory

for i in ${input_files[@]}
do
    echo "Calculating coverage at Chromosom ends (50kb) for file: ${i}"
    name=$(basename ${i})

    echo "Calculating coverage along chrom. Q arms:"
    bedtools multicov -bams ${i} -bed ${chrEnds_Q} > ${outdir}/${name%.bam}_chrQ_multicov_genmap600bp.bed
    
    echo "Calculating coverage along chrom. P arms:"
    bedtools multicov -bams ${i} -bed ${chrEnds_P} > ${outdir}/${name%.bam}_chrP_multicovgenmap600bp.bed
 
done

echo "Done."

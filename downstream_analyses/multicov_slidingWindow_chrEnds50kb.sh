#!/bin/bash -l

sliding_win_chrP="/CHM13/slidingWindows_bedtools_chromEndCoordinates_50kb_chrParms_w100_s50.bed"
sliding_win_chrQ="/CHM13/slidingWindows_bedtools_chromEndCoordinates_50kb_chrQarms_w100_s50.bed"

inputdir=$(realpath $1)
chrEnds_P="/CHM13/CHM13_hg38chrY_chromEndCoordinates_50kb_chrParms.genome"
chrEnds_Q="/CHM13/CHM13_hg38chrY_chromEndCoordinates_50kb_chrQarms.genome"
outdir=$(realpath $2)

input_files=$(ls ${inputdir}/*.bam | tr "\n" " ")

for i in ${input_files[@]}
do
    echo "Calculating coverage at Chromosome ends (50kb) for file: ${i}"
    name=$(basename ${i})

    echo "Calculating coverage along chrom. Q arms:"
    samtools view -b ${bam} -L ${chrEnds_P} | genomeCoverageBed -ibam stdin -bg | bedtools map -b stdin -a ${sliding_win_chrQ} -c 4 -o count

    bedtools multicov -bams ${outdir}/${name%.bam}_rev.bam -bed ${chrEnds_Q} > ${outdir}/${name%.bam}_rev_chrQ_multicov.bed
    
    echo "Calculating coverage along chrom. P arms:"
    bedtools multicov -bams ${outdir}/${name%.bam}_rev.bam -bed ${chrEnds_P} > ${outdir}/${name%.bam}_rev_chrP_multicov.bed

    bedtools multicov -bams ${outdir}/${name%.bam}_fwd.bam -bed ${chrEnds_Q} > ${outdir}/${name%.bam}_fwd_chrQ_multicov.bed
 
done

echo "Done."

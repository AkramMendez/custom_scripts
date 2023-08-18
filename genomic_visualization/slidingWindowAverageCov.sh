#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 03:00:00
#SBATCH -J slidingWindowCov


module load bioinfo-tools
module load BEDTools

# Calculate mean coverage of sliding windows for IP and input samples.
# Input files: sorted sliding window map from reference genome previously generated with Bedtools.
# i.e: bedtools makewindows -g ${chromsizes} -w 100 -s 50 > bedtools.sliding.windows.20size.10step.hybGenome.hg38.masked.humanSTF500.ecoli.bed
# bedgraph file= bedgraph file for calculating the average from.


sliding_window_bed="bedtools.sliding.windows.20size.10step.hybGenome.hg38.masked.humanSTF500.ecoli.bed"

bedgraph=$1

echo "Calculating mean coverage over sliding windows:"


# From sliding window method from Wenqui et al 2021 "METTL3 regulates heterochromatin in mouse embryonic stem cells":
# m6 A peak calling was performed using a ‘sliding window’ method slightly modified from a previous study[22] . 

#In brief, read numbers of IP and input were calculated on every 25-bp window across the genome. 

bedtools multicov -bams $(ls | grep -P "COVID_UK.*.withdup.uniq.sorted.fwd.bam$" | tr "\n" " ") -bed ../../reference_genomes/wuhCor1/wuhCor1.slidingWin.width25.step10.bed > bedtools.multicov.covidUK.withdup.uniq.m6a.input.fwd

bedtools map -a ${sliding_window_bed} -b <(bedtools bamtobed -i ${bam}) -c 4 -o count -null 0 > ${bam%.bam}.genomecov.slidwin25.bedgraph

bedtools intersect -wao -a ${bed_ip} -b ${bed_input} | awk 'BEGIN{FS=OFS="\t"}{print $1,$2,$3,$4/$8}' > ratio.IP.input.bedgraph

#Windows with normalized IP/input density ≥ 2 and Fisher’s exact test P < 0.05 were selected.
use bedgraphTobigwig to convert mapped coverage to bigwig ant then 
use bigwigCompare to get the ratio files from IP and input in bedgraph output and the use bedtools fisher to compare both.


#Adjacent windows were merged using bedtools (v2.29.2)41 and merged regions with size ≥100 bp were determined as m6A peaks.

echo "Done."

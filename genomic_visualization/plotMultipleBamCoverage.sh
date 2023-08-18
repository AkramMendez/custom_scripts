#!/bin/bash -l
# Script for generating a coverage track in bigwig format from multiple bamfiles BAM file using deepTools bamCoverage function.

scaling_factors_file=$(realpath $1) #Scaling factors to normalize the coverage tracks, the file must specify the BAM path and its corresponding scaling factor in a tab-separated format
n=$2 # Number of processors

plotBamCoverage(){

	echo "plotting bamCoverage for sample ${bam}"
	bamCoverage -b ${bam} -o bigwigs/${bam%.bam}.rev.bw --scaleFactor ${sf} --binSize 1 -p ${n}
}

while IFS=$'\t' read -r bam sf 
do
	plotBamCoverage ${bam} ${sf} ${n}

done < <(cat ${scaling_factors_file} | tail -n  +2)
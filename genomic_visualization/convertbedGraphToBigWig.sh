


module load bioinfo-tools
module load ucsc-utilities

bedgraph=$1
chrom_sizes=$2
outdir=$3
outbed=${bedgraph/.bdg/.sortbed.bdg}
bigwig=${bedgraph/.sortbed.bdg/.bw}

sort -k1,1 -k2,2n ${bedgraph} > ${outbed}
#bedSort ${bedgraph} ${outbed}

bedGraphToBigWig ${outbed} ${chrom_sizes} ${bigwig}
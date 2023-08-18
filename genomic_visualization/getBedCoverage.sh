#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 06:00:00
#SBATCH -J bedCoverage

module load bioinfo-tools
module load deepTools

IP=$1
input=$2
outdir=$3
out_IP=${IP%.bam}.bdg
out_input=${input%.bam}.bdg

n=$(echo "$(($SLURM_NTASKS-1))")

echo "Generating coverage bedgraph for IP and input samples:"
bamCoverage --bam ${IP} --outFileName ${outdir}/${out_IP} --outFileFormat bedgraph --normalizeUsing RPKM --binSize 5 --numberOfProcessors ${n} --verbose
bamCoverage --bam ${input} --outFileName ${outdir}/${out_input} --outFileFormat bedgraph --normalizeUsing RPKM --binSize 5 --numberOfProcessors ${n} --verbose

echo "Generating log2 ratio bedgraph IP vs input"
bamCompare --bamfile1 ${IP} --bamfile2 ${input} --outFileName ${outdir}/${out_IP}_vs_${out_input}.log2ratio --outFileFormat bedgraph --normalizeUsing RPKM --binSize 5 --numberOfProcessors ${n} --verbose --scaleFactorsMethod None

echo "Done."



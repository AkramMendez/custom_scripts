#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 18:00:00
#SBATCH -J metpeakCsh

module load R_packages/4.0

gtf="hybrid.hg38.gencondev36.chr.subtel.ecoli.gtf"

outdir="/metpeakResultsV2/"
metpeak_script="metpeakAnalysis.R"

#cd ${TMPDIR}

IP=$1
input=$2
name=$3
#gtf=$4
#outdir=$4

Rscript --vanilla --verbose ${metpeak_script} ${IP} ${input} ${name} ${gtf} ${outdir}

#
#echo "Resulting files appear as:"
#ls -l ${TMPDIR}
#
#echo "Moving results to ${outdir}"
#
#mv ${TMPDIR}/${name} ${outdir}
#
echo "Done."


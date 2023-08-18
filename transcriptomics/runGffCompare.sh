#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 03:00:00
#SBATCH -J gffcompare


module load bioinfo-tools

#This script compares GFF or GTF files produced by Stringtie assembly against a reference genome GTF and generates several statistics to evaluate the quality of the assembly.

gtf_ref=$(realpath $1)
gtf_stringtie=$(realpath $2)
outdir=$(realpath $3)

gffcompare="/domus/h1/amendez/private/gffcompare/gffcompare"

cd ${outdir}

echo "Comparing ${gtf_stringtie} to ${gtf_ref}"

${gffcompare} -R -r ${gtf_ref} -o strtcmp ${gtf_stringtie}

echo "Done"

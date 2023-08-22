#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 16:00:00
#SBATCH -J fetchSRA


module load bioinfo-tools
module load sratools
module load gnuparallel

runids=$1
outdir=$2

echo "Fetching data from SRA accession list and spliting FASTQ files"

cat ${runids} | parallel --eta --verbose  "fastq-dump -O ${outdir} --split-files -F {}"

echo "Done"
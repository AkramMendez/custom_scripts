#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 24:00:00
#SBATCH -J Stringtie2

#This scripts performs a de novo transcriptome assembly from long-reads using StringTie2
#https://github.com/gpertea/stringtie

module load bioinfo-tools

stringtie="/domus/h1/amendez/private/stringtie/stringtie"

inputdir=$(realpath $1) # Input directory with BAM files
outdir=$(realpath $2) # Output directory

cd ${SNIC_TMP}

echo "Copying BAM files to ${SNIC_TMP}"

cp ${inputdir}/*.bam .
cp ${inputdir}/*.bai .

for bam in $(ls *.bam)
do
    echo "Performing de novo transcriptome assembly with StringTie2 for ${bam}"
    outfile=$(echo "${bam%.bam}_txm_Stringtie2.2.1.gtf")
    ${stringtie} -L -o ${outfile} ${bam} -A -v
done


echo "Copying assembled files to ${outdir}"
cp *.gtf ${outdir}
cp *.tab ${outdir}

echo "Done."
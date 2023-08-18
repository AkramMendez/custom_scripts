#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 8
#SBATCH -t 04:00:00
#SBATCH -J annotPeaks


module load bioinfo-tools
module load HOMER

ref="CHM13_T2T.fa"
gff3="CHM13_T2T.gff3"

inputdir=$(realpath $1)
outdir=$(realpath $2)

cd ${SNIC_TMP}

cp ${inputdir}/*.narrowPeak

for peak in $(ls *.narrowPeak)
do

    name=$(basename ${peak} | sed 's/.narrowPeak//g')

    echo "Annotating peaks for ${name}"

    annotatePeaks.pl <(awk 'BEGIN{FS=OFS="\t"}{if($4~/fwd/){gsub(".","+",$6); print}}' ${peak}) ${ref} -strand + -annStats ${outdir}/annStats_${name}_fwd.txt -gff3 ${gff3} -gid -cpu ${SLURM_NTASKS} > ${outdir}/annotatePeaks_${name}_fwd.txt

    annotatePeaks.pl <(awk 'BEGIN{FS=OFS="\t"}{if($4~/rev/){gsub(".","-",$6); print}}' ${peak}) ${ref} -strand - -annStats ${outdir}/annStats_${name}_rev.txt -gff3 ${gff3} -gid -cpu ${SLURM_NTASKS} > ${outdir}/annotatePeaks_${name}_rev.txt

done
echo "Done."

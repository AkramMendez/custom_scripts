#!/bin/bash -l
#SBATCH -A naiss2023-22-153
#SBATCH -p core -n 8
#SBATCH -t 12:00:00
#SBATCH -J annotPeaks


module load bioinfo-tools
module load HOMER

ref="/crex/proj/nb_project/private/genomes/grch38/GRCh38.primary_assembly.genome.fa"
gff3="/proj/nb_project/private/genomes/grch38/gencode.v38.annotation.gtf"


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
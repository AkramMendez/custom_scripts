#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 10
#SBATCH -t 12:00:00
#SBATCH -J trimmingReadsBGI

module load bioinfo-tools
module load TrimGalore/0.6.1

inputdir="/crex/proj/nb_storage/private/ncc_project/raw_data/clean_bgi_tncc_be2"
outdir="/crex/proj/nb_storage/private/ncc_project/trimmed_reads"


#After sequencing quality control, the first 10 and the last five bases were trimmed

cd ${inputdir}

echo "Trimming fastq files"
for i in ${names[@]}; do R1=$(echo ${i}"_1.fq.gz"); R2=$(echo ${i}"_2.fq.gz"); echo "Processing paired samples: ${R1} and ${R2}" ; trim_galore ${i}/${R1} ${i}/${R2} --paired --length 20 --cores ${SLURM_NTASKS} -o ${outdir} ;done

echo "Done."

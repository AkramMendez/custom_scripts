#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 10
#SBATCH -t 08:00:00
#SBATCH -J disambiguate

# Activate the xengsort environment (I installed NGS disambiguate in this env.) https://github.com/AstraZeneca-NGS/disambiguate
conda activate xengsort
module load bioinfo-tools

inputdir_mouse=$(realpath $1)
inputdir_human=$(realpath $2)
outdir=$(realpath $3)

cd ${SNIC_TMP}

mkdir mouse
mkdir human

echo "Copying mouse alignments"
cp ${inputdir_mouse}/*.bam ./mouse
cp ${inputdir_mouse}/*.bai ./mouse

echo "Copying human alignments"

cp ${inputdir_human}/*.bam ./human
cp ${inputdir_human}/*.bai ./human

samples_mouse=$(ls ./mouse/*.bam)

samples_human=$(ls ./human/*.bam)

sample_ids=$(ls ./mouse/*.bam | sed -E 's/(mouse|.bam)//g' | sort | uniq | tr "\n" " ")

for sample in ${sample_ids[@]}
do
    echo "Running disambiguate on sample: ${sample} mouse and human alignments"

    ngs_disambiguate -s ${sample} -o ${outdir} -a star ./human/${sample}_human.bam ./mouse/${sample}_mouse.bam

done

echo "Done."


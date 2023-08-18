#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 10
#SBATCH -t 24:00:00
#SBATCH -J xengsort

conda activate xengsort
module load bioinfo-tools

mouse=$(realpath $1)
human=$(realpath $2)
inputdir=$(realpath $3)
outdir=$(realpath $4)

cd ${SNIC_TMP}

echo "Copying samples to ${SNIC_TMP}"

cp ${inputdir}/*.fq.gz .

echo "Building xengsort index"

xengsort index  -H <(zcat ${mouse}) -G <(zcat ${human}) -n 4_500_000_000 -T ${SLURM_NTASKS}  myindex.h5

samples=$(ls *.fq.gz | sed -E 's/_(1|2).fq.gz//g' | sort | uniq | tr "\n" " ")

for i in ${samples[@]}
do
    echo "Classifying reads for ${i} condition "
    xengsort classify --index myindex.h5  --fastq <(zcat ${i}_1.fq.gz)  --pairs <(zcat ${i}_2.fq.gz) -T ${SLURM_NTASKS} --out ${outdir}/${i}_xengsort_results

done

echo "Moving files to ${outdir}"

cp *.h5 ${outdir}
cp *.fq ${outdir}

echo "Done."


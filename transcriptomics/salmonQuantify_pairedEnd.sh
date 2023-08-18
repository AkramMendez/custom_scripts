#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 10
#SBATCH -t 12:00:00
#SBATCH -J salmonQuant

module load bioinfo-tools
module load Salmon/1.9.0

gtf=$(realpath $1)
salmon_idx=$(realpath $2)
input_reads=$(realpath $3)
outdir=$(realpath $4)

cd ${outdir}

for i in $(find ${input_reads} -mindepth 1 -type d); do name=$(basename ${i}); echo "--- ${name} ----"; sample=$(basename ${i}/*.gz | sed -E 's/_(1|2).fq.gz//g'); echo "Processing sample: ${sample}"; R1=${i}/${sample}_1.fq.gz; R2=${i}/${sample}_2.fq.gz; echo "Reads $R1 and $R2"; salmon quant -i ${salmon_idx} -l A -g ${gtf} --validateMappings --gcBias --minScoreFraction 0.5 -1 ${R1} -2 ${R2} -o salmon_quant_${name} ; done

echo "Done."

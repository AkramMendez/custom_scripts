#!/bin/bash -l
#SBATCH -A naiss2023-22-153
#SBATCH -p core -n 8
#SBATCH -t 08:00:00
#SBATCH -J salmonQuant

module load bioinfo-tools
module load Salmon/1.4.0

#Generate "gentrome" (genome + transcriptome):
#cat gencode.v38.transcripts.fa GRCh38.primary_assembly.genome.fa > gentrome.salmon.GRCh38.fa

#Make decoys file, containing the chromosome names to map the transcriptome sequences to:
#grep "^>" GRCh38.primary_assembly.genome.fa | cut -d " " -f 1 > decoys.txt
#sed -i.bak -e 's/>//g' decoys.txt

inputdir=$(realpath $1)
outdir=$(realpath $2)

salmon_idx="/crex/proj/nb_project/private/genomes/grch38/salmon_index_grch38"
gtf="/crex/proj/nb_project/private/genomes/grch38/gencode.v38.annotation.gtf"

cd ${SNIC_TMP}

echo "Copying files to ${SNIC_TMP}"
cp ${inputdir}/*.gz .

for read in $(ls *.gz)
do
    name=$(basename ${read} | sed 's/_trimmed.fq.gz//g')

    echo "Quantifying sample ${name}"
    salmon quant -i ${salmon_idx} -l A -g ${gtf} --validateMappings -r ${read} -o salmon_quant_${name} -p ${SLURM_NTASKS}

done

echo "Copying output file to ${outdir}"

cp -R salmon_quant_* ${outdir}

echo "Done."
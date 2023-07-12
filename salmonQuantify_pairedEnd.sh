#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 02:00:00
#SBATCH -J salmonQuant

module load bioinfo-tools
module load Salmon/1.4.0

# Script for quantifying paired-end reads using Salmon with GRCh38 human reference genome

#Generate "gentrome" (genome + transcriptome):
#cat gencode.v38.transcripts.fa GRCh38.primary_assembly.genome.fa > gentrome.salmon.GRCh38.fa

#Make decoys file, containing the chromosome names to map the transcriptome sequences to:
#grep "^>" GRCh38.primary_assembly.genome.fa | cut -d " " -f 1 > decoys.txt
#sed -i.bak -e 's/>//g' decoys.txt

salmon_idx=$(realpath $1) # Path to the salmon index
gtf=$(realpath $2) # Path to GTF annotation file (i.e gencode.v38.annotation.gtf)
read1=$(realpath $3) # Path to read 1 fastq.gz file
read2=$(realpath $4) # Path to read 2 fastq.gz file
name=$5 # output sample name

echo "Quantifying reads ${read1} and ${read2} using Salmon"

salmon quant -i ${salmon_idx} -l A -g ${gtf} --validateMappings -1 ${read1} -2 ${read2} -o salmon_quant_${name} -p ${SLURM_NTASKS}

echo "Done."
#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 24:00:00
#SBATCH -J minimap2.22

module load bioinfo-tools

minimap2="/domus/h1/amendez/private/minimap2/minimap2"

inputdir=$(realpath $1)
outdir=$(realpath $2)
logdir=$(realpath $3)
kmer_size=5
ref_fasta="/ref/CHM13/chm13.draft_v1.1.transcriptome.gffread.fasta"

# Script for mapping reads to the T2T human reference genome using Minimap v2.22, optimized for mapping reads into highly-repetitive regions
# The k-mer size parameter can be set to small values for increased sensitivity (https://github.com/lh3/minimap2)

mapReadsMinimap2(){
	ref=$1
	fastq=$2
	kmer_size=$3
	samfile=$4

	${minimap2} -ax splice -uf -k${kmer_size} -t ${SLURM_NTASKS} --secondary=no ${ref_fasta} ${fastq} > ${samfile} 2>> ${logdir}
}

echo "Mapping reads to ${ref_fasta} reference using Minimap v2.22:"
mapReadsMinimap2 ${ref_fasta} ${fastq} ${kmer_size} ${samfile}

echo "Done."
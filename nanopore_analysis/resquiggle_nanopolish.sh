#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p node -n 32
#SBATCH -t 12:00:00
#SBATCH -J nanopolish

module load bioinfo-tools
module load nanopolish

fastq=$(realpath $1) #Path to basecalled fastQ reads
fast5=$(realpath $2) #Path to fast5 files
bam=$(realpath $3)  #Path to BAM alignment
genome=$(realpath $4)   #Path to reference genome fasta file
summary=$(realpath $5)  #Path to sequencing_summary.txt file
outdir=$(realpath $6)   #Path to output directory
name=$7 #Name for saving the output eventalign file.

# Script for aligning Nanopore events to a reference using Nanopolish (https://nanopolish.readthedocs.io/en/latest/quickstart_eventalign.html)
# The eventaligned can the be used for downstream analysis such as m6A detection using the xPore pipeline (https://xpore.readthedocs.io/en/latest/).

cd $SNIC_TMP

echo "Index with Nanopolish"
nanopolish index -d ${fast5} ${fastq}

echo "Nanopolish eventalign"
nanopolish eventalign --reads ${fastq} \
--bam ${bam} \
--genome ${genome} \
--signal-index \
--scale-events \
--summary ${summary} \
--progress \
--threads ${SLURM_NTASKS} >> ${outdir}/${name}_eventalign.txt

echo "Copying index to ${outdir}"
cp *.index ${outdir}
cp *.fai ${outdir}
cp *.gzi ${outdir}
cp *.readdb ${outdir}

echo "Done."

#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 24:00:00
#SBATCH -J minimap2.22

module load bioinfo-tools
module load samtools

minimap2="/domus/h1/amendez/private/minimap2/minimap2"

fastq=$(realpath $1)
name=$2
outdir=$(realpath $3)
kmer_size=14

#Align to CHM13 transcriptome
ref_fasta="/ref/CHM13/chm13.draft_v1.1.transcriptome.gffread.fasta"

#Align to CHM13 reference genome:
#ref_fasta="/ref/CHM13/chm13.draft_v1.1.fasta"

cd ${SNIC_TMP}

#${minimap2} -ax splice -uf -k${kmer_size} -t ${SLURM_NTASKS} --secondary=no ${ref_fasta} ${fastq} > ${name}.chm13.txm.k14.sam
echo "Mapping with Minimap v2.22:"

${minimap2} -ax map-ont -uf -L -k${kmer_size} -t ${SLURM_NTASKS} --secondary=no ${ref_fasta} ${fastq} > ${name}.chm13.txm.k14.sam

samtools view -Sb ${name}.chm13.txm.k14.sam | samtools sort -o ${name}.chm13.txm.k14.sorted.bam -

samtools index ${name}.chm13.txm.k14.sorted.bam

echo "Copying BAM and index files to ${outdir}"

cp *.bam ${outdir}

cp *.bai ${outdir}

echo "Done."


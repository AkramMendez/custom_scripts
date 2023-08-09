#!/bin/bash -l
#SBATCH -A snic2022-22-85
#SBATCH -p core -n 12
#SBATCH -t 2-00:00:00
#SBATCH -J minimap2.24

# Script for mapping basecalled Nanopore long-reads to a reference transcriptome or genome using Minimap v 2.24

module load bioinfo-tools
module load samtools
module load deepTools/3.3.2

conda activate /proj/nb_project/private/conda_envs/minimap2/

fastq=$(realpath $1) # Path to basecalled FastQ file
name=$2 # Sample name
outdir=$(realpath $3) # Output directory
kmer_size=$4 # k-mer size
gx_or_tx=$5 # Map to genome or transcriptome ["gx", "tx"]
secondary=$6 # Keep secondary mappings ["yes", "no"]

cd ${SNIC_TMP}

if [[ ${gx_or_tx} == "tx" ]] && [[ ${secondary} == "yes" ]]
then
    ref_fasta="chm13.v1.1.transcriptome.gffread.fasta"
    
    echo "Align to CHM13 reference transcriptome ${ref_fasta}"
    echo "Mapping with Minimap v2.24:"
    
    minimap2 -ax map-ont -uf \
    -k${kmer_size} \
    -t ${SLURM_NTASKS} \
    --secondary=yes ${ref_fasta} ${fastq} > ${name}.chm13.txm.${kmer_size}.sam

    samtools view -bS -F 20 ${name}.chm13.txm.${kmer_size}.sam | samtools sort -o ${name}.chm13.txm.${kmer_size}.mapped.plus.bam -@ ${SLURM_NTASKS}

    echo "Indexing BAM file"
    samtools index ${name}.chm13.txm.${kmer_size}.mapped.plus.bam -@ ${SLURM_NTASKS}

    echo "Making bigwig track for ${name}.chm13.txm.${kmer_size}.mapped.plus.bam"
    bamCoverage -b ${name}.chm13.txm.${kmer_size}.mapped.plus.bam \
    -p ${SLURM_NTASKS} \
    --outFileFormat bigwig \
    --binSize 10 \
    -o ${name}.chm13.txm.${kmer_size}.mapped.plus.bs10.bw

else
    ref_fasta="chm13.v1.1.fasta"
    
    echo "Align to CHM13 reference genome ${ref_fasta}"
    echo "Mapping with Minimap v2.24:"

    minimap2 -ax map-ont -uf \
    -k${kmer_size} \
    -t ${SLURM_NTASKS} \
    --secondary=yes ${ref_fasta} ${fastq} > ${name}.chm13.genome.${kmer_size}.sam
    
    samtools view -bS -F 20 ${name}.chm13.genome.${kmer_size}.sam | samtools sort -o ${name}.chm13.genome.${kmer_size}.mapped.plus.bam -@ ${SLURM_NTASKS}

    echo "Indexing BAM file"
    samtools index ${name}.chm13.genome.${kmer_size}.mapped.plus.bam -@ ${SLURM_NTASKS}

    echo "Making bigwig track for ${name}.chm13.genome.${kmer_size}"
    bamCoverage -b ${name}.chm13.genome.${kmer_size}.mapped.plus.bam \
    -p ${SLURM_NTASKS} \
    --outFileFormat bigwig \
    --binSize 10 \
    -o ${name}.chm13.genome.${kmer_size}.mapped.plus.bs10.bw
fi


echo "Copying BAM and index files to ${outdir}"

cp *.bam ${outdir}
cp *.bw ${outdir}
cp *.bai ${outdir}

echo "Done."
#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 06:00:00
#SBATCH -J bowtie_aln

module load bioinfo-tools
module load gcc/10.2.0 star/2.7.2b


STAR --runThreadN 10 \
--runMode genomeGenerate \
--genomeDir /proj/snic2020-15-304/private/star_index_hg38 \
--genomeFastaFiles /domus/h1/amendez/genomes/gencode/grch38/GRCh38.primary_assembly.genome.fa \
--sjdbGTFfile /domus/h1/amendez/genomes/gencode/grch38/gencode.v36.chr_patch_hapl_scaff.annotation.gtf \
--sjdbOverhang 100

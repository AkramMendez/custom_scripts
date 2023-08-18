#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 10
#SBATCH -t 08:00:00
#SBATCH -J sortmerna

module load bioinfo-tools
module load SortMeRNA/3.0.3

reads=$1
outdir=$2

cd ${SNIC_TMP}

for fastq in $(ls ${reads}/*.fq.gz)
do

	name=$(basename ${fastq})

	sortmerna --ref $SORTMERNA_DBS/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5.8s-database-id98:$SORTMERNA_DBS/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5s-database-id98:$SORTMERNA_DBS/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNA_DBS/index/silva-euk-18s-id95:$SORTMERNA_DBS/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMERNA_DBS/index/silva-euk-28s-id98 \
	--reads-gz ${fastq} \
	--aligned ./rRNA_reads_${name%.fq.gz}.fastq \
	--other ./non_rRNA_reads_${name%.fq.gz}.fastq \
	--fastx --log -v -a 2
	
done

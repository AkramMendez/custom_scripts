#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 06:00:00
#SBATCH -J filter_rRNA




module load SortMeRNA

The rRNA databases provided with SortMeRNA have been indexed.

    rfam-5.8s-database-id98.fasta
    rfam-5s-database-id98.fasta
    silva-arc-16s-id95.fasta
    silva-arc-23s-id98.fasta
    silva-bac-16s-id90.fasta
    silva-bac-23s-id98.fasta
    silva-euk-18s-id95.fasta
    silva-euk-28s-id98.fasta

The Fasta files may be found at $SORTMERNA_DBS/rRNA_databases/
and the indices may be found at $SORTMERNA_DBS/index/, with the
same name as the corresponding Fasta file minus the '.fasta' suffix.
For example:

sortmerna --ref $SORTMERNA_DBS/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5.8s-database-id98: \ 
$SORTMERNA_DBS/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5s-database-id98: \ 
$SORTMERNA_DBS/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNA_DBS/index/silva-euk-18s-id95: \
$SORTMERNA_DBS/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMERNA_DBS/index/silva-euk-28s-id98 \
--reads ${reads} \
--fastx \
--aligned ${prefix}
--log -v


sortmerna --ref $SORTMERNA_DBS/rRNA_databases/rfam-5.8s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5.8s-database-id98:$SORTMERNA_DBS/rRNA_databases/rfam-5s-database-id98.fasta,$SORTMERNA_DBS/index/rfam-5s-database-id98:$SORTMERNA_DBS/rRNA_databases/silva-euk-18s-id95.fasta,$SORTMERNA_DBS/index/silva-euk-18s-id95:$SORTMERNA_DBS/rRNA_databases/silva-euk-28s-id98.fasta,$SORTMERNA_DBS/index/silva-euk-28s-id98 --reads ${reads} --fastx --aligned -a 4 --log -v

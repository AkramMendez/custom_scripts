#!/bin/bash -l
#SBATCH -A snic2020-15-304
#SBATCH -p core -n 8
#SBATCH -t 01:00:00
#SBATCH -J annotatePeaksHomer

reference="hybridHg38SubtelEcoli_v3/hybGenome.hg38.masked.humanSTF500.ecoli.fa"
gtf="hybGenome.hg38.masked.humanSTF500.ecoli.1based.gtf"

peaks=$1
sample=$2

annotatePeaks.pl ${peaks} ${reference} -gtf ${gtf} > annotatePeaks.${sample}.${peaks%.bed}.txt

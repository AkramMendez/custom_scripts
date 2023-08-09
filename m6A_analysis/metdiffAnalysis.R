#!/usr/bin/env Rscript

library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(IRanges)
library(tidyverse)
library(MetDiff, lib.loc="/domus/h1/amendez/R/x86_64-pc-linux-gnu-library/4.0")


args = commandArgs(trailingOnly=TRUE)


gtf <- args[1]
outdir <- args[2]

setwd(outdir)

knockdown.rev<-"knockdown.hisat2.nodup.sorted.negrev.bam"
knockdown.input.rev<-"knockdown.input.hisat2.nodup.sorted.negrev.bam"
control.rev<-"control.hisat2.nodup.sorted.negrev.bam"
control.input.rev<-"control.input.hisat2.nodup.sorted.negrev.bam"

metdiff(GENE_ANNO_GTF = gtf ,IP_BAM = control.rev,INPUT_BAM = control.input.rev, TREATED_IP_BAM = knockdown.rev,TREATED_INPUT_BAM = knockdown.input.rev, EXPERIMENT_NAME = "KD_vs_Control")






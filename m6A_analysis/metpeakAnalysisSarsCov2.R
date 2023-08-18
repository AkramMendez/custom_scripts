#!/usr/bin/env Rscript

library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(IRanges)
library(tidyverse)
library(MeTPeak,lib.loc="/domus/h1/amendez/R/x86_64-pc-linux-gnu-library/4.0")



args = commandArgs(trailingOnly=TRUE)

ip <- args[1]
input <- args[2]
name<- args[3]
gtf<- args[4]
outdir<-args[5]

setwd(outdir)

metpeak(GENE_ANNO_GTF = gtf ,IP_BAM = ip,INPUT_BAM = input, EXPERIMENT_NAME = name)
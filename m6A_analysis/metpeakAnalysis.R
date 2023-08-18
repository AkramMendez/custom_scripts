#!/usr/bin/env Rscript

library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(IRanges)
library(tidyverse)
library(MeTPeak,lib.loc="/domus/h1/amendez/R/x86_64-pc-linux-gnu-library/4.0")


args = commandArgs(trailingOnly=TRUE)

IP<-args[1]
input<-args[2]
name<-args[3]
gtf=args[4]
workdir=args[5]

setwd(workdir)

cat("Working dir:",getwd())

message("Satrting MeTPeak analysis")

metpeak(GENE_ANNO_GTF = gtf ,IP_BAM = IP,INPUT_BAM = input, EXPERIMENT_NAME = name)

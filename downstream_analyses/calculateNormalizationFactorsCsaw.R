#!/usr/bin/env Rscript

library(csaw)
library(tidyverse)
library(data.table)

args = commandArgs(trailingOnly=TRUE)

bamdir=normalizePath(args[1])

bamfiles<-list.files(bamdir,pattern = ".bam$")

binnedCounts<-windowCounts(bam.files = bamfiles,width = 1000,ext = 75)

norm_factors<-normFactors(binnedCounts)


names(norm_factors$norm.factors)<-gsub("_S\\d+_R\\d+_\\d+.*.bam$","",norm_factors$bam.files,perl = T)
norm_factors$norm.factors

cat("TMM-Normalization factors:\n",norm_factors)
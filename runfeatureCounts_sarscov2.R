#!/usr/bin/env Rscript

library(Rsubread)


args = commandArgs(trailingOnly=TRUE)

inputdir <- args[1]
design_table_scov2 <- args[2]
threads <- args[3]


bamfiles<-list.files("/mappings/covid19_viral_m6a_input",pattern = "markdup.bam$")
bamfiles

sarscov2_annotation<-""



featureCounts_scov2<-featureCounts(bamfiles, annot.ext=sarscov2_annotation, isGTFAnnotationFile = TRUE, nthreads = threads)

counts_scov2<-featureCounts_scov2$counts

DESeqDataSetFromMatrix
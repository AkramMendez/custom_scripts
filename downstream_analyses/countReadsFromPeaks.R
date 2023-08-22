library(Rsubread)
library(tidyverse)


# This script calculates raw reads located within coordinates of interest (such as peak coordinates) in a BAM file. 
# The counts are obtained using featureCounts and coordinates provided in SAF format.
# Raw counts can be further normalized to CPM, TPM, RPKM, etc. for downstream analysis.

path <- "/bamfiles_directory/"

bamfiles <- list.files(path, pattern = "*.bam$", full.names=TRUE)

colnamesPeaks <- c("chr","start","end","name","score","strand","fold_enrichment","-log10pvalue","-log10qvalue","peak.point.source")
peaks <- fread("/peaks/peaks.narrowPeak",col.names = colnamesPeaks)

# Convert narrowPeaks to SAF format
peaks_SAF <- peaks %>% dplyr::select(GeneID=name,Chr=chr,Start=start, End=end,Strand=strand) %>% mutate(Strand=if_else(grepl("_rev_",GeneID),"-","+"))


countsFromPeaks <- featureCounts(files = bamfiles,
                                 annot.ext = peaks_SAF,
                                 isGTFAnnotationFile = FALSE, 
                                 useMetaFeatures = TRUE, 
                                 isPairedEnd = FALSE,
                                 minOverlap = 1, 
                                 primaryOnly = FALSE, 
                                 strandSpecific = 1, 
                                 verbose = TRUE)

countsFromPeaks$counts %>% head()
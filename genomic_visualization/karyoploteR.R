library(tidyverse)
library(data.table)
library(karyoploteR)
library(GenomicRanges)
library(rtracklayer)
library(regioneR)
library(plyranges)


fai <- fread("/CHM13/CHM13_hg38chrY_ecoliK12.fa.fai", col.names=c("chr", "chrlen", "x","y","z"))
fai <- fai %>% filter(!grepl("U00096.2|chrM",chr))
# Data frame with chromosome lengths
chroms <- data.table(fai$chr,0,fai$chrlen, gieStain="geng")
cytobands <- fread("cytoBandT2TCHM13_ucsc.txt", col.names=c("chr","start","end","name","gieStain"))
chm13.genome <- toGRanges(chroms)

chrEnds <- import.bed("CHM13_hg38chrY_chromEndCoordinates_50kb_chrPandQarms.genome")

kd_fwd <- subsetByOverlaps(import.bw("/bigwig/merged_normSpikeIn_2/kd.mapped.nodup.spikeInNorm.bs10.ratio.fwd.bw"),chrEnds)

kd_rev <- subsetByOverlaps(import.bw("/bigwig/merged_normSpikeIn_2/kd.mapped.nodup.spikeInNorm.bs10.ratio.rev.bw"),chrEnds)

names(mcols(chrQ30kb.fwd.bw)) <- "value"

colNames <-c("chr","start","end","name","score","strand","fold_enrichment","-log10pvalue","-log10qvalue","summit_position")

#kd chrP file size is 0 because no peaks found for chrP in kd
kd_peaksChrP <- fread("bedtools_intersect_KD_noTreat_ext75_peaks_all_50kb_chrParms.txt",sep = "\t", header = F)[,1:10]


kd_peaksChrQ<-fread("bedtools_intersect_KD_noTreat_ext75_peaks_all_50kb_chrQarms.txt", sep = "\t", header = F)[,1:10]
colnames(kd_peaksChrQ) <- colNames

kd_chrQP <- mergeRegions(kd_peaksChrP,kd_peaksChrQ)

control_peaksChrP <- fread("bedtools_intersect_Control_noTreat_noTreat_ext75_peaks_all_50kb_chrParms.txt", sep = "\t", header = F)[,1:10]

control_peaksChrQ <- fread("bedtools_intersect_Control_noTreat_noTreat_ext75_peaks_all_50kb_chrQarms.txt", sep = "\t", header = F)[,1:10]

colnames(control_peaksChrP) <- colNames
colnames(control_peaksChrQ) <- colNames

control_peaksChrP <- makeGRangesFromDataFrame(control_peaksChrP,keep.extra.columns = T)

control_peaksChrQ <- makeGRangesFromDataFrame(control_peaksChrQ,keep.extra.columns = T)

control_chrQP <- plyranges::bind_ranges(control_peaksChrP,control_peaksChrQ)

#kp <- plotKaryotype(genome=chm13.genome,cytobands = cytobands, ideogram.plotter = NULL,chromosomes = "chr1",zoom = toGRanges("chr1:1-50000"))
#zoom = toGRanges("chr1:1-2000000") chromosomes = "chr1"
#Cytobands interesect with chrom ends 50kb
cytobands_chrEnd50kb <- plyranges::join_overlap_inner(chrEnds,makeGRangesFromDataFrame(cytobands,keep.extra.columns = T))

#xpore datasets

chrmsWithRipPeaks <- unique(
  c(
    control_chrQP %>% as.data.frame() %>% dplyr::select(seqnames) %>% distinct() %>% unlist(),
    kd_chrQP %>% as.data.frame() %>% dplyr::select(seqnames) %>% distinct() %>% unlist()
  )
)

modRatesAtChrEnds_control_drach <- modRatesAtChrEnds_control %>% as.data.frame() %>% filter(mod_rate_CTRLSH.REP1 > 0 & grepl("..AC.",kmer,perl = T)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)

modRatesAtChrEnds_kd_drach <- modRatesAtChrEnds_kd %>% as.data.frame() %>% filter(mod_rate_KD.REP1 > 0 & grepl("..AC.",kmer,perl = T)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)

#log10 scaled genome coordinates
#chm13.genome <- chm13.genome %>% as.data.frame() %>% mutate(start=log10(start +1), end=log10(end +1))
#cytobands <- cytobands %>% as.data.frame() %>% mutate(start=log10(start +1), end=log10(end +1))

pdf("karyotypePlot_peaks.pdf")

kp <- plotKaryotype(genome=chm13.genome %>% filter(seqnames %in% chrmsWithRipPeaks),cytobands = cytobands %>% filter(chr %in% chrmsWithRipPeaks), ideogram.plotter = NULL,plot.type = 1,main = "\n m6A-RIP peaks \n (Control-sh and KD)",plot.params = )
#"\n xPore NNACN and MeRIP peaks (Ctrl-sh)"
kpAddCytobands(kp,lwd=0.33,clipping = FALSE)

#kpPlotRegions(kp,modRatesAtChrEnds_control_drach, col = "blue",num.layers = 2,avoid.overlapping = F)
kpPlotRegions(kp,control_chrQP, col = "red",num.layers = 2,avoid.overlapping = F)
kpPlotRegions(kp,kd_chrQP, col = "blue",num.layers = 2,avoid.overlapping = F)

dev.off()


kp2 <- plotKaryotype(genome=chm13.genome,cytobands = cytobands, ideogram.plotter = NULL,plot.type = 1,main = "\n xPore NNACN and MeRIP peaks (KD)")

kpAddCytobands(kp2,lwd=0.33)

kpPlotRegions(kp2,modRatesAtChrEnds_kd, col = "green",num.layers = 2,avoid.overlapping = F)
kpPlotRegions(kp2,modRatesAtChrEnds_control %>% filter(mod_rate_CTRLSH-REP1 > 0), col = "orange",num.layers = 2,avoid.overlapping = F)

kpPlotBigWig(kp, data="/bigwig/chrQ_reads/fwd/kd_mapped.nodup.sorted.bam.spikeInNorm.chrQ.fwd.bs10.bw",ymax = 0.5,r0 = 0.1,r1=0.5,col = "steelblue",clipping = F,)

kpPlotBigWig(kp, data="/bigwig/chrQ_reads/fwd/kd_mapped.nodup.sorted.bam.spikeInNorm.chrQ.fwd.bs10.bw",ymax = 1,r0 = 0.1,r1=1,col = "steelblue")


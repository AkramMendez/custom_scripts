library(data.table)
library(tidyverse)

# Script for processing m6Anet results obtained from Nanopore long-read data to extract the most enriched predicted modified k-mers in the whole CHM13 T2T human reference transcriptome.
# Read the transcript coordinates from a BED file
tx_coords <- fread("transcript_coordinates_chm13.draft_v1.1.gene_annotation.v4.gffread.bed", sep = "\t", col.names = c("chr","start","end","tx_id","strand"))

# Read the control and treatment m6aNet results:
m6anet_ctrl <- fread("m6anetResults_ctrl_t2tTxm.csv", sep=",", header = T)

m6anet_treat <- fread("m6anetResults_treat_t2tTxm.csv", sep=",", header = T)

# Filter m6ANet results
m6anet_ctrl_filt <- m6anet_ctrl %>% filter(n_reads > 5, probability_modified > 0, mod_ratio > 0.25)

m6anet_treat_filt <- m6anet_treat %>% filter(n_reads > 5, probability_modified > 0.9, mod_ratio > 0.5)

 
sites_ctrl <- left_join(m6anet_ctrl_filt, tx_coords, by=c("transcript_id"="tx_id")) %>% 
  mutate(startSite=start+transcript_position,endSite=start+transcript_position+6) %>%
  dplyr::select(chr,startSite,endSite,kmer,mod_ratio,strand,start,end)


sites_treat <- left_join(m6anet_treat_filt, tx_coords, by=c("transcript_id"="tx_id")) %>% 
  mutate(startSite=start+transcript_position,endSite=start+transcript_position+6) %>% 
  dplyr::select(chr,startSite,endSite,kmer,mod_ratio,strand)

fwrite(sites_ctrl,"../m6aNet_results/sites_m6aNet_ProbAbove0.5_ModRateAbove0.5_ctrl.bed", col.names = F, sep="\t")

fwrite(sites_treat,"../m6aNet_results/sites_m6aNet_ProbAbove0.5_ModRateAbove0.5_treat.bed",col.names = F, sep="\t")

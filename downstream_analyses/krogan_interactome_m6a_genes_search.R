library(tidyverse)
library(openxlsx)


res_all_lfcShrink<-fread("diffexp_lfcShrink.tsv",header = T,sep="\t")

cov2genes<-c("S","ORF8","ORF7a","ORF7b","ORF6","ORF3a","ORF1ab","ORF10","N","M","E","ORF1a")

degs_df<-res_all_lfcShrink %>% filter(!(gene %in% cov2genes) & !grepl("gene-",gene) & abs(log2FoldChange))

highConfInt_cov2<- read.xlsx("/public_data/interactome_MS_cov2_krogan2020/SuppTable2_high_confidence_interactors.xlsx", colNames = T)
allInteractions_cov2<-read.xlsx("/public_data/interactome_MS_cov2_krogan2020/Supplementary_Table_1_scoring_all_baits_and_proteins.xlsx", colNames = T)


highConfInt_N <- highConfInt_cov2 %>% filter(Bait=="SARS-CoV2 N") %>% dplyr::select(PreyGene) %>% distinct() %>% unlist()

allInteractions_N <- allInteractions_cov2 %>% filter(Bait=="SARS-CoV2 N") %>% dplyr::select(PreyGene) %>% distinct() %>% unlist()

#Are they differentially expressed?
#268 interactors from 583 N protein interactors are differentially expressed.
diffexp_Interactors_N <- res_all_lfcShrink %>% filter(gene %in% allInteractions_N & abs(log2FoldChange) >1 & padj < 0.01 & !is.na(padj))

deg_genes_Interactors_N <- res_all_lfcShrink %>% filter(gene %in% allInteractions_N & abs(log2FoldChange) >1 & padj < 0.01 & !is.na(padj)) %>% dplyr::select(gene) %>% distinct() %>% unlist()

#How many have DEU?
#59 are DEG and DEU
dex_filt %>% filter(gene %in% allInteractions_N) %>% dim()

res_all_lfcShrink %>% filter(gene %in% deg_genes_Interactors_N & log2FoldChange > 1 & padj < 0.01 & !is.na(padj)) %>% dim()
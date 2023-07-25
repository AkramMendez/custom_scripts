




# Reported genes after infection took from the supplementary data of Wyler et al., 2021 
#(Transcriptomic profiling of SARS-CoV-2 infected human cell lines identifies HSP90 as target for COVID-19 therapy)
# https://ars.els-cdn.com/content/image/1-s2.0-S258900422100119X-mmc1.pdf

#The authors compared upregulated genes after SARS-CoV-2 infection in Calu-3 and Caco-2 cells with other public data from incetion with Cilliated nasopharyngeal swabs (Chua et al., 2020),
#bronchoalveolar lavages (BALs) (Liao et al., 2020),and transcriptome studies of normal human bronchial epithelial (NHBE) cells, 
#A549 cells with and without ACE2 expression, and also Calu-3 cells upon infection with SARS-CoV-2 (Blanco-Melo et al., 2020) 


# Genes exhibiting changes both in patient epithelial cells and Calu-3 infected cells.
sharedCaluAndEpithelialCellsPatients <- c("IFIT3","IRF7","IFIT1","NFKBIA","RSAD2","IFIT2","ATF3","MX2","OAS1","XAF1","CXCL10","MX1","PPP1R15A","IFI44L","IFITM1","JUN","CCL2","BST2","FOSB","KLF6","IFITM3","MT2A","IFITM2","ISG15","IFI27","IFI6","DDIT3","TXNIP","ARRDC3")

# Upregulated genes both in Caco-2 and Calu-3 after infection
signifChangedCacoAndCalu <- c("ARRDC3","ATF3","CCL2","CXCL11","CXCL10","DUSP1","EGR1","FOS","IFIT1","IFIT2","IFIT3","IFIT5","IFNL1","IFNB1","IL6","ISG15","JUN","JUNB","NFKBIA","OAS1","OAS3","PPP1R15A","TNF","TTR","TXNIP")


# Genes upregulated in patient cilliated nasopharyngeal swabs and Calu-3 12 hpi
sharedCilliatedAndCalu <- c("IFITM1","CCL5","ICAM1","FOS","IFITM3","PSMB9","LTB","IFITM2","PLAUR","PTGS2","CXCL1","ISG20","IL32","PNRC1","CD83","ARRDC3","TXNIP","NFKBIA","TNFAIP3","PPP1R15A")

# Upregulated in patient Secretory nasopharyngeal swabs and Calu-3 12 hpi
sharedSecretoryAndCalu <- c("CXCL3","IFITM1","BST2","IFITM3","PSMB9","CXCL1","ISG20","IFI27","IFI6","ARRDC3","TXNIP","PPP1R15A")

# Upregulated in A549 series 2 and 5 (Blano-Melo) and Calu-3 12 hpi
sharedA549series2AndCalu <- c("ATF3","CCL2","CMPK2","CXCL10","CXCL11","CYP1A1","IDO1","IFI44L","IFIT1","IFIT2","IFIT3","IFNB1","IFNL1","IFNL2","IFNL3","IL6","MX1","MX2","RSAD2","TNF","TNFAIP3","TRIM22","XAF1","ARRDC3","PPP1R15A","TXNIP")

sharedA549series5AndCalu <- c("ATF3","CCL2","CMPK2","CXCL10","CXCL11","CYP1A1","IDO1","IFIT1","IFIT2","IFIT3","IFNB1","IFNL1","IFNL2","IFNL3","IL6","MX2","RSAD2","TNF","TNFAIP3","TRIM22","TXNIP","ARRDC3","PPP1R15A")


intersect(sharedCaluAndEpithelialCellsPatients,diffexpUp_wu)
intersect(sharedCilliatedAndCalu,diffexpUp_wu)
intersect(sharedSecretoryAndCalu,diffexpUp_wu)
intersect(sharedA549series2AndCalu,diffexpUp_wu)

intersect(sharedCaluAndEpithelialCellsPatients,diffexpDown_wu)

intersect(sharedCaluAndPatients,diffexpUp_uk)
intersect(sharedCaluAndPatients,diffexpUp_sa)

sharedAcrossPatientAndCellLines <- unique(c(sharedCaluAndEpithelialCellsPatients,sharedCilliatedAndCalu,sharedSecretoryAndCalu,sharedA549series2AndCalu,sharedA549series5AndCalu))



res_all_lfcShrink %>%
  dplyr::select(gene,log2FoldChange,padj,contrast) %>% distinct() %>% 
  filter(abs(log2FoldChange)>1 & padj < 0.01 &! is.na(padj) & gene %in% sharedCaluAndEpithelialCellsPatients) %>% 
  mutate(comparison="Shared genes Epithelial Patient Cells and Calu-3") %>% 
  ggplot(aes(gene,log2FoldChange, fill=contrast)) + 
  geom_bar(stat = "identity", position = position_dodge2(width = 0.9,preserve = "single")) + 
  labs(title = "Shared genes in Patien Epithelial Cells and Calu-3 after infection") + 
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_viridis_d(begin = 0.5)


res_all_lfcShrink %>%
  dplyr::select(gene,log2FoldChange,padj,contrast) %>% distinct() %>% 
  filter(abs(log2FoldChange)>1 & padj < 0.01 &! is.na(padj) & gene %in% sharedSecretoryAndCalu) %>% 
  mutate(comparison="Shared genes Secretory Patient Cells and Calu-3") %>% 
  ggplot(aes(gene,log2FoldChange, fill=contrast)) + 
  geom_bar(stat = "identity", position = position_dodge2(width = 0.9,preserve = "single")) + 
  labs(title = "Shared genes in Patien Secretory Cells and Calu-3 after infection") + 
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_viridis_d(begin = 0.5)


res_all_lfcShrink %>%
  dplyr::select(gene,log2FoldChange,padj,contrast) %>% distinct() %>% 
  filter(abs(log2FoldChange)>1 & padj < 0.01 &! is.na(padj) & gene %in% sharedAcrossPatientAndCellLines) %>% 
  mutate(comparison="Shared genes in Patients and Cell lines") %>% 
  ggplot(aes(gene,log2FoldChange, fill=contrast)) + 
  geom_bar(stat = "identity", position = position_dodge2(width = 0.9,preserve = "single")) + 
  labs(title = "Shared genes in Patients and Cell lines") + 
  theme(axis.text.x = element_text(angle = 45)) +
  scale_fill_viridis_d(begin = 0.5)

# Differentially expressed genes after Wu, UK or SA infection
res_all_lfcShrink %>%
  dplyr::select(gene,log2FoldChange,padj,contrast) %>% distinct() %>% 
  filter(abs(log2FoldChange)>1 & padj < 0.01 &! is.na(padj) & gene %in% sharedAcrossPatientAndCellLines) %>% dplyr::select(gene) %>% distinct() %>% unlist() -> tmp


#Differentially expressed genes from Blanco-Melo et al,.2020:
# "Differentially expressed genes (DEGs) were characterized for each sample (|L2FC| > 1, p-adjusted-value < 0.05)"
meloData <- readxl::read_xlsx("~/MondalLab/BEA21P074_Roshan_Vaid/plots/cov2_figures_2022/differentially_expressed_genes_patients_covid19_Blanco-Melo2020.xlsx", sheet = 2, col_types = c("chr","num","num","num","num","num","chr"))

meloData <- meloData %>% filter(status=="OK")
meloData <- meloData %>% filter(abs(log2FoldChange) >2 & padj < 0.05 &! is.na(padj))

meloData %>% filter(log2FoldChange > 2) %>% dplyr::select(Gene_name) %>% distinct() %>% unlist() -> upregMelo

meloData %>% filter(log2FoldChange < -2) %>% dplyr::select(Gene_name) %>% distinct() %>% unlist() -> downregMelo


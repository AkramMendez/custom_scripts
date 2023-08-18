library(tidyverse)
library(openxlsx)
library(kableExtra)

cov2genes<-c("S","ORF8","ORF7a","ORF7b","ORF6","ORF3a","ORF1ab","ORF10","N","M","E","ORF1a")

#VIRMA=KIAA1429, HAKAI=CBLL1, "HNRNPG"=RBMX
#m6A-related proteins: https://jhoonline.biomedcentral.com/articles/10.1186/s13045-020-00872-8/tables/1

writers<-sort(c("METTL3","METTL14","METTL16","KIAA1429","CBLL1","WTAP","RBM15","RBM15B","ZC3H13","SPEN"))
readers<-sort(c("YTHDC1","YTHDC2","YTHDF1","YTHDF2","YTHDF3","HNRNPA2B1","HNRNPC","RBMX","IGF2BP1","IGF2BP2","IGF2BP3"))
erasers<-sort(c("FTO","ALKBH5"))
stress_granule<-c("ATXN2","ATXN2L","BICD1","CIRBP","D1PAS1","DDX3X","DDX6","DYNC1H1","G3BP1","G3BP2","GRB7","OGFOD1","PQBP1","PRKAA1","PRKAA2","PRRC2C","PUM2","RPS23","STYXL1","UBAP2L")

degs_df<-res_all_lfcShrink %>% filter(!(gene %in% cov2genes) & !grepl("gene-",gene) & abs(log2FoldChange))

highConfInt_cov2<- read.xlsx("interactome_MS_cov2_krogan2020/SuppTable2_high_confidence_interactors.xlsx", colNames = T)
allInteractions_cov2<-read.xlsx("interactome_MS_cov2_krogan2020/Supplementary_Table_1_scoring_all_baits_and_proteins.xlsx", colNames = T)


highConfInt_N <- highConfInt_cov2 %>% filter(Bait=="SARS-CoV2 N") %>% dplyr::select(PreyGene) %>% distinct() %>% unlist()

allInteractions_N <- allInteractions_cov2 %>% filter(Bait=="SARS-CoV2 N") %>% dplyr::select(PreyGene) %>% distinct() %>% unlist()

rnaInteractome <- read.xlsx("interactome_cov2_munschauer2021/SCoV2-RNA-protein-Atlas-Supplemental-Tables/Supplementary_Table_1.xlsx", colNames = T)

### From Schmidt et al. 
# "We identified 276 proteins with a positive log2 fold change. Of these, 57 were significantly
# enriched (adjusted P < 0.05, purifications, two-tailed t-test), which we subsequently defined as the set of core
# SARS-CoV-2 RNA interacting proteins (Fig. 1c). Additionally, we also
# defined an expanded SARS-CoV-2 RNA interactome using a relaxed false discovery rate (FDR) of less than 20%
# (Fig. 1c)."
###

#276 set log2 TMT ratio (SARS-CoV-2/RMRP) > 0
rnaInteractome <- rnaInteractome %>% filter(logFC.SCoV2.over.RMRP >0)

# Core SARS-CoV-2 RNA interactome: log2 TMT ratio (SARS-CoV-2/RMRP) > 0 and adjusted p-value < 0.05
rnaInteractome %>% filter(logFC.SCoV2.over.RMRP >0 & adj.P.Val.SCoV2.over.RMRP < 0.05)

# Expanded RNA interactome:log2 TMT ratio (SARS-CoV-2/RMRP) > 0 and adjusted p-value < 0.2
#104 host + SARS-CoV-2 proteins set, 118 proteins
rnaInteractome <- rnaInteractome %>% filter(logFC.SCoV2.over.RMRP >0 & adj.P.Val.SCoV2.over.RMRP < 0.2)


rnaInteractome <- rnaInteractome %>% mutate(geneSymbol=if_else(is.na(geneSymbol),id,geneSymbol))

rnaInteractome$geneSymbol %>% as.data.frame() %>% fwrite(.,"list_geneNames_rnaInteractome_munschauer2021.txt", sep="\t",col.names = F)

#STRING-db interactions m6A readers, no more than 50 interactors first shell:
stringInteractions_m6areaders <- fread("string_interactions_experiment_noMore500_physical_m6Areaders_interacting_SARSCov2_RNA_interactome_Schmidt2021.tsv", sep="\t", header=T)

#Nodes on string-db network
length(unique(union(stringInteractions_m6areaders$node2,stringInteractions_m6areaders$node1)))

string_nodes <- unique(union(stringInteractions_m6areaders$node2,stringInteractions_m6areaders$node1))

#Note: from the reported 276 interactors, two have duplicated assigned gene symbols, I pasted an additional digit to distiguish those duplicated ids.
pdf("vennDiagram_m6areaders_STRING_vs_SARS-COV2_RNAInteractomeSchmidt2021.pdf")

ggvenn(list("SARS-CoV-2 RNA-prot"=c(unique(rnaInteractome$geneSymbol),"YBX3.1","CNBP.1"),"String-db"=string_nodes), text_size = 8,show_percentage = F,stroke_color = NA, auto_scale=TRUE)

dev.off()

#Common interactors m6A readers and SARS-CoV2-RNA interactome proteins
intersect(unique(rnaInteractome$geneSymbol),unique(c(stringInteractions_m6areaders$node1,stringInteractions_m6areaders$node2)))

# m6a-related proteins in RNA interactome ---------------------------------

rnaInteractome %>% filter(geneSymbol %in% c(readers,writers,erasers)) %>% dplyr::select(geneSymbol)


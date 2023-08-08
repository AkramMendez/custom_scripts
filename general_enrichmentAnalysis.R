require(DESeq2)
require(tidyverse)
require(pheatmap)
require(data.table)
require(ggtext)
require(factoextra)
require(pheatmap)
require(clusterProfiler)
require(enrichR)
require(AnnotationDbi)
require(org.Hs.eg.db)
require(tidyverse)
require(enrichplot)
require(msigdbr)

##############################################################################################################################
### Perform GSEA analysis
##############################################################################################################################
# Preparing MsigDB datasets -----------------------------------------------
#### See https://www.gsea-msigdb.org/gsea/msigdb/human/genesets.jsp
# **c2.all        -curated gene sets
# **c2.cpg        -chemical and genetic pertubations
# **c2.cp         -canonical pathways
# **c2.wikipath   -wikipathways
# **c5.all        -ontology gene sets: BP, CC, MF, Human Phenotype Ontology
# **c6.all        -oncogenic signature sets
# **c8.all        -cell type signature gene sets
# **h.hallmarks   -hallmark gene sets
##############################################################################################################################

# HALLMARK
#Hallmark gene sets (50 total gene sets)
msigdbH <- msigdbr(species = "Homo sapiens",category = "H") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()
# C2
#Curated gene sets (6290 gene sets, this can be split into the following subcategories, see below:)
    #msigdbC2<- msigdbr(species = "Homo sapiens",category = "C2") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()

    #Curated gene sets WikiPathways (615 gene sets)
    msigdbWP<- msigdbr(species = "Homo sapiens",category = "C2", subcategory = "CP:WIKIPATHWAYS") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()
    #chemical and genetic pertubations (3368 gene sets)
    msigdbCGP<- msigdbr(species = "Homo sapiens",category = "C2", subcategory = "CGP") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()
    #canonical pathways (29 gene sets)
    msigdbCP<- msigdbr(species = "Homo sapiens",category = "C2", subcategory = "CP") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()
    # Reactome (1604 gene sets)
    msigdbREACTOME<- msigdbr(species = "Homo sapiens",category = "C2", subcategory = "CP:REACTOME") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()
    # PID (Pathway Interaction Database) (196 gene sets)
    msigdbPID<- msigdbr(species = "Homo sapiens",category = "C2", subcategory = "CP:PID") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()
    # KEGG (186 gene sets)
    msigdbKEGG<- msigdbr(species = "Homo sapiens",category = "C2", subcategory = "CP:KEGG") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()
    # BIOCARTA (292 gene sets)
    msigdbBIOCARTA <- msigdbr(species = "Homo sapiens",category = "C2", subcategory = "CP:BIOCARTA") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()
    
#Regulatory target gene sets (3731 gene sets)
#msigdbC3<- msigdbr(species = "Homo sapiens",category = "C3") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()

# Computational gene sets (858 gene sets)
# msigdbC4<- msigdbr(species = "Homo sapiens",category = "C4") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()

# Gene ontology gene sets:
# "GO:BP" "GO:CC" "GO:MF" "HPO"
# GO Biologucal process (7481)
msigdbC5_GOBP <- msigdbr(species = "Homo sapiens",category = "C5", subcategory = "GO:BP") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()
# GO MF (1708 gene sets)
msigdbC5_GOMF <- msigdbr(species = "Homo sapiens",category = "C5", subcategory = "GO:MF") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()
# GO CC (996 gene sets)
msigdbC5_GOCC <- msigdbr(species = "Homo sapiens",category = "C5", subcategory = "GO:CC") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()
# HPO (Human Phenotype Ontology)
msigdbC5_HPO <- msigdbr(species = "Homo sapiens",category = "C5", subcategory = "HPO") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()

#################
#Convert named ranked list to Entrez gene symbols, the vector needs to have the ranked values and has to be named using Entrez IDs  
ranked_list_entrez <- ranked_list

names(ranked_list_entrez) <- mapIds(org.Hs.eg.db,keys=names(ranked_list), column="SYMBOL",keytype="ENTREZID", multiVals="first")


# Run GSEA analysis -------------------------------------------------------------
#Run this in Uppmax:

gseaH <-        GSEA(sort(ranked_list_symbols, decreasing=TRUE), TERM2GENE = msigdbH,  pvalueCutoff = 0.05)
gseaWP  <-      GSEA(sort(ranked_list_symbols, decreasing=TRUE), TERM2GENE = msigdbWP,  pvalueCutoff = 0.05)
gseaCGP <-      GSEA(sort(ranked_list_symbols, decreasing=TRUE), TERM2GENE = msigdbCGP,  pvalueCutoff = 0.05)
gseaCP <-       GSEA(sort(ranked_list_symbols, decreasing=TRUE), TERM2GENE = msigdbCP,  pvalueCutoff = 0.05)
gseaREACTOME <- GSEA(sort(ranked_list_symbols, decreasing=TRUE), TERM2GENE = msigdbREACTOME,  pvalueCutoff = 0.05)
gseaPID <-      GSEA(sort(ranked_list_symbols, decreasing=TRUE), TERM2GENE = msigdbH,  pvalueCutoff = 0.05)
gseaKEGG <-     GSEA(sort(ranked_list_symbols, decreasing=TRUE), TERM2GENE = msigdbH,  pvalueCutoff = 0.05)
gseaBIOCARTA <- GSEA(sort(ranked_list_symbols, decreasing=TRUE), TERM2GENE = msigdbH,  pvalueCutoff = 0.05)
gseaC5_GOBP <-  GSEA(sort(ranked_list_symbols, decreasing=TRUE), TERM2GENE = msigdbH,  pvalueCutoff = 0.05)
gseaC5_GOMF <-  GSEA(sort(ranked_list_symbols, decreasing=TRUE), TERM2GENE = msigdbH,  pvalueCutoff = 0.05)
gseaC5_GOCC <-  GSEA(sort(ranked_list_symbols, decreasing=TRUE), TERM2GENE = msigdbH,  pvalueCutoff = 0.05)
gseaC5_HPO <-   GSEA(sort(ranked_list_symbols, decreasing=TRUE), TERM2GENE = msigdbH,  pvalueCutoff = 0.05)


#Oncogenic signature gene sets (189 gene sets)
msigdbC6 <- msigdbr(species = "Homo sapiens",category = "C6") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()
#Cell type signature gene sets (671 gene sets)
msigdbC8 <- msigdbr(species = "Homo sapiens",category = "C8") %>% dplyr::select(gs_name,human_gene_symbol) %>% distinct()





degsLFC_kd1 <- diffGenesSalmon_CellType1_kd1 %>% arrange(desc(log2FoldChange)) %>% dplyr::select(log2FoldChange) %>% unlist()
names(degsLFC_kd1) <- diffGenesSalmon_CellType1_kd1 %>% arrange(desc(log2FoldChange)) %>% dplyr::select(gene) %>% unlist()

degsLFC_kd2 <- diffGenesSalmon_CellType1_kd2 %>% arrange(desc(log2FoldChange)) %>% dplyr::select(log2FoldChange) %>% unlist()
names(degsLFC_kd2) <- diffGenesSalmon_CellType1_kd2 %>% arrange(desc(log2FoldChange)) %>% dplyr::select(gene) %>% unlist()

degsLFC_kd1_entrez <- mapIds(org.Hs.eg.db,keys=names(degsLFC_kd1),column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")


# Enrichment Diseases -----------------------------------------------------

enrichNCG_CellType1_kd1 <- enrichNCG(degsLFC_kd1_entrez, pvalueCutoff = 0.05, readable = T)
enrichDisgenet_CellType1_kd1 <- enrichDGN(degsLFC_kd1_entrez, pvalueCutoff = 0.05, readable = T)

enrichDO_CellType1_kd1 <- enrichDO(degsLFC_kd1_entrez,pvalueCutoff = 0.05, readable = T)

enrichDO_CellType1_kd1@result %>% dplyr::filter(grepl("neuroblastoma",Description))

enrcichDisgenet_CellType1_kd1@result %>% dplyr::filter(grepl("Neurob",Description) & pvalue < 0.05)


# Load gsea GO analysis run in Uppmax:
lfcskd1 <- lfcShrinSalmon_CellType1_kd1 %>% dplyr::filter(!is.na(log2FoldChange)) %>% arrange(log2FoldChange) %>% dplyr::select(log2FoldChange) %>% unlist()
names(lfcskd1) <- lfcShrinSalmon_CellType1_kd1 %>% dplyr::filter(!is.na(log2FoldChange)) %>% arrange(log2FoldChange) %>% dplyr::select(gene) %>% unlist()

lfcskd1_entrez <- mapIds(org.Hs.eg.db, names(lfcskd1),column="ENTREZID",keytype="SYMBOL",multiVals="first")


names(lfcskd1) <- lfcskd1_entrez

lfcskd1 <- lfcskd1[which(!is.na(names(lfcskd1)))]

lfcskd1 <- lfcskd1[which(duplicated(names(lfcskd1))==FALSE)]

#lfcskd1 <- lfcskd1[rank(sort(lfcskd1, decreasing = T),ties.method = "max")]

save(lfcskd1,file="/tables/lfcskd1.Rda")

lfcskd2 <- lfcShrinSalmon_CellType1_kd2 %>% dplyr::filter(!is.na(log2FoldChange)) %>% arrange(log2FoldChange) %>% dplyr::select(log2FoldChange) %>% unlist()
names(lfcskd2) <- lfcShrinSalmon_CellType1_kd2 %>% dplyr::filter(!is.na(log2FoldChange)) %>% arrange(log2FoldChange) %>% dplyr::select(gene) %>% unlist()

lfcskd2_entrez <- mapIds(org.Hs.eg.db, names(lfcskd2),column="ENTREZID",keytype="SYMBOL",multiVals="first")


names(lfcskd2) <- lfcskd2_entrez

lfcskd2 <- lfcskd2[which(!is.na(names(lfcskd2)))]

lfcskd2 <- lfcskd2[which(duplicated(names(lfcskd2))==FALSE)]


save(lfcskd2,file="/tables/lfcskd2.Rda")


# GSEA GO processes analysis run in Uppmax ---------------------------------------------



#gseBP_CellType1kd1 <- gseGO(sort(lfcskd1,decreasing=T), org.Hs.eg.db, ont="BP", pvalueCutoff=0.05, eps=0,seed=123)
#gseBP_CellType1kd1 <- setReadable(gseBP_CellType1kd1,org.Hs.eg.db)

#gseMF_CellType1kd1 <- gseGO(sort(lfcskd1,decreasing=T), org.Hs.eg.db, ont="MF", pvalueCutoff=0.05, eps=0,seed=123)
#gseMF_CellType1kd1 <- setReadable(gseMF_CellType1kd1,org.Hs.eg.db)

#gseCC_CellType1kd1 <- gseGO(sort(lfcskd1,decreasing=T), org.Hs.eg.db, ont="CC", pvalueCutoff=0.05, eps=0,seed=123)
#gseCC_CellType1kd1 <- setReadable(gseCC_CellType1kd1,org.Hs.eg.db)

# gseBP_CellType1kd2 <- gseGO(sort(lfcskd2,decreasing=T), org.Hs.eg.db, ont="BP", pvalueCutoff=0.05, eps=0,seed=123)
# gseBP_CellType1kd2 <- setReadable(gseBP_CellType1kd2,org.Hs.eg.db)
# 
# gseMF_CellType1kd2 <- gseGO(sort(lfcskd2,decreasing=T), org.Hs.eg.db, ont="MF", pvalueCutoff=0.05, eps=0,seed=123)
# gseMF_CellType1kd2 <- setReadable(gseMF_CellType1kd2,org.Hs.eg.db)
# 
# gseCC_CellType1kd2 <- gseGO(sort(lfcskd2,decreasing=T), org.Hs.eg.db, ont="CC", pvalueCutoff=0.05, eps=0,seed=123)
# gseCC_CellType1kd2 <- setReadable(gseCC_CellType1kd2,org.Hs.eg.db)


load("/tables/gseBP_CellType1kd1.Rda")
load("/tables/gseBP_CellType1kd2.Rda")
load("/tables/gseMF_CellType1kd1.Rda")


#load("/tables/gseMF_CellType1kd2.Rda")
#load("/tables/gseCC_CellType1kd1.Rda")
#load("/tables/gseCC_CellType1kd2.Rda")

# Enrichr analysis kd1-sh -------------------------------------------------


enrichrWP_CellType1_kd1 <- enrichr(names(degsLFC_kd1),databases = "WikiPathway_2021_Human")
enrichrMsigH_CellType1_kd1 <- enrichr(names(degsLFC_kd1),databases = "MSigDB_Hallmark_2020")
enrichrMsigOnco_CellType1_kd1 <- enrichr(names(degsLFC_kd1),databases = "MSigDB_Oncogenic_Signatures")
enrichrReactome_CellType1_kd1 <- enrichr(names(degsLFC_kd1),databases = "Reactome_2016")
enrichrKEGG_CellType1_kd1 <- enrichr(names(degsLFC_kd1),databases = "KEGG_2021_Human")
#neural crest differentiation: c("MSX2","OLIG3","PAX3","ZIC5","CDH7","CDH6","BMP4","FGF8","GBX2","FGF19","HEY2","SNAI2","HES1","MSX1","FGFR3","HES5")

#gseaGOBP2_DEGs_CellType1_kd1 <- gseGO(geneList = degsLFC_kd1,OrgDb =org.Hs.eg.db,keyType = "SYMBOL",ont="BP",pvalueCutoff = 0.05)

mydotplot<-function(df,title, showCategory,x){
  Count<-x
  df@result$GeneRatio<-sapply(df@result$GeneRatio,function(x){eval(parse(text=x))})
  df %>% slot("result") %>% arrange(pvalue) %>% filter(pvalue<0.05) %>% head(showCategory) %>% ggplot(aes(reorder(str_wrap(Description,70),Count),Count,color=pvalue,size=GeneRatio)) + 
    geom_point() + coord_flip() + scale_color_gradient(low="red", high="blue") + 
    theme_light()+ 
    theme(panel.background = element_blank(), panel.grid = element_blank(), aspect.ratio = 1.5) +
    labs(title=title) + guides(color = guide_colorbar(reverse = TRUE)) + xlab("") + ylab("Count") -> p
  return(p)
}


# GO enrichment kd1-sh ----------------------------------------------------
goBP_CellType1_kd1<- enrichGO(gene = names(degsLFC_kd1),OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% simplify() %>% gofilter(.,level = c(7:10))

goMF_CellType1_kd1<- enrichGO(gene = names(degsLFC_kd1),OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% simplify() %>% gofilter(.,level = c(7:10))

goCC_CellType1_kd1<- enrichGO(gene = names(degsLFC_kd1),OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% simplify() %>% gofilter(.,level = c(7:10))


reactomePA_CellType1_kd1 <- enrichPathway(degsLFC_kd1_entrez,pvalueCutoff = 0.05, organism = "human", readable = T)


# Enrichment plots kd1 and kd2 --------------------------------------------

pdf("/plots/enrichmentAnalysis_NCCday6_kd1sh.pdf", width = 12, height = 14)

mydotplot(goBP_CellType1_kd1 , showCategory=30, x="Count", title="GO:Biological Process \n DEGs CellType1 kd1-sh")
mydotplot(goMF_CellType1_kd1 , showCategory=30, x="Count", title="GO:Molecular function \n DEGs CellType1 kd1-sh")
mydotplot(goCC_CellType1_kd1 , showCategory=30, x="Count", title="GO:Cellular component \n DEGs CellType1 kd1-sh")

mydotplot(reactomePA_CellType1_kd1 , showCategory=30, x="Count", title="GO:Cellular component \n DEGs CellType1 kd1-sh")
mydotplot(enrichDO_CellType1_kd1, showCategory = 30,x="Count", title = "Enrichment Disease Ontology \n DEGs CellType1 kd1-sh")
mydotplot(enrichNCG_CellType1_kd1, showCategory = 30,x="Count", title = "Enrichment NCG Network of Cancer Genes \n DEGs CellType1 kd1-sh")
mydotplot(enrichDisgenet_CellType1_kd1, showCategory = 30,x="Count", title = "Enrichment DisGeNet \n DEGs CellType1 kd1-sh")

plotEnrich(enrichrWP_CellType1_kd1$WikiPathway_2021_Human %>% filter(P.value < 0.05),numChar = 70, title = "Enrichr: WikiPathways \n DEGs CellType1 kd1-sh")
plotEnrich(enrichrKEGG_CellType1_kd1$KEGG_2021_Human %>% filter(P.value < 0.05),numChar = 70, title = "Enrichr: KEGG 2021 \n DEGs CellType1 kd1-sh")
plotEnrich(enrichrReactome_CellType1_kd1$Reactome_2016 %>% filter(P.value < 0.05),numChar = 70, title = "Enrichr: Reactome \n DEGs CellType1 kd1-sh")
plotEnrich(enrichrMsigH_CellType1_kd1$MSigDB_Hallmark_2020 %>% filter(P.value < 0.05),numChar = 70, title = "Enrichr: MSigDB Hallmark gene sets \n DEGs CellType1 kd1-sh")
plotEnrich(enrichrMsigOnco_CellType1_kd1$MSigDB_Oncogenic_Signatures %>% filter(P.value < 0.05),numChar = 70, title = "Enrichr: MSigDB Oncogenic Signatures \n DEGs CellType1 kd1-sh")

gseaplot2(gseaH_CellType1_kd1,geneSetID =gseaH_CellType1_kd1@result$ID, title = paste("GSEA MSigDB [H-Hallmark gene sets] \n CellType1 kd1-sh \n", gseaH_CellType1_kd1@result$ID))
gseaplot2(gseaC2_CellType1_kd1,geneSetID =gseaC2_CellType1_kd1@result$ID[1:5], title = "GSEA MSigDB [C2-Curated gene sets] \n CellType1 kd1-sh")
gseaplot2(gseaC3_CellType1_kd1,geneSetID =gseaC3_CellType1_kd1@result$ID[1:5], title = "GSEA MSigDB [C3-Regulatory target gene sets] \n CellType1 kd1-sh")
gseaplot2(gseaC4_CellType1_kd1,geneSetID =gseaC4_CellType1_kd1@result$ID[1:5], title = "GSEA MSigDB [C4-Computational gene sets] \n CellType1 kd1-sh")
gseaplot2(gseaC8_CellType1_kd1,geneSetID =gseaC8_CellType1_kd1@result$ID[1:5], title = "GSEA MSigDB [C8-Cell type signature gene sets] \n CellType1 kd1-sh")

gseaplot2(gseBP_kd1 %>% simplify(), geneSetID = 1:5,title = "GSEA GO:Biological Process \n CellType1 kd1-sh")
gseaplot2(gseMF_kd1, geneSetID = gseMF_kd1@result$ID,title="GSEA GO:Molecular function \n CellType1 kd1-sh")
gseaplot2(gseCC_kd1, geneSetID = gseCC_kd1@result$ID, title= paste("GSEA GO:Cellular component \n CellType1 kd1-sh \n", gseCC_kd1@result$ID))

dev.off()

enrichTables_kd1 <- list(
  "GO_BP"=goBP_CellType1_kd1@result %>% filter(pvalue < 0.05),
  "GO_MF"=goMF_CellType1_kd1@result %>% filter(pvalue < 0.05),
  "GO_CC"=goCC_CellType1_kd1@result %>% filter(pvalue < 0.05),
  "Reactome"=reactomePA_CellType1_kd1@result %>% filter(pvalue < 0.05),
  "DiseaseOnt"=enrichDO_CellType1_kd1@result %>% filter(pvalue < 0.05),
  "NetCancerGenes"=enrichNCG_CellType1_kd1@result %>% filter(pvalue < 0.05),
  "DisGenNet"=enrichDisgenet_CellType1_kd1@result %>% filter(pvalue < 0.05),
  "WikiPaths"=enrichrWP_CellType1_kd1$WikiPathway_2021_Human %>% filter(P.value < 0.05),
  "KEGG"=enrichrKEGG_CellType1_kd1$KEGG_2021_Human %>% filter(P.value < 0.05),
  "Enrichr_Reactome"=enrichrReactome_CellType1_kd1$Reactome_2016 %>% filter(P.value < 0.05),
  "Enrichr_HallmarksMsigdb"=enrichrMsigH_CellType1_kd1$MSigDB_Hallmark_2020 %>% filter(P.value < 0.05),
  "Enrichr_OncoMsigdb"=enrichrMsigOnco_CellType1_kd1$MSigDB_Oncogenic_Signatures %>% filter(P.value < 0.05),
  "GSEA_Hallmark_MSigdb"=gseaH_CellType1_kd1@result %>% filter(pvalue < 0.05),
  "GSEA_Curated_MSigdb"=gseaC2_CellType1_kd1@result %>% filter(pvalue < 0.05),
  "GSEA_RegTarget_MSigdb"=gseaC3_CellType1_kd1@result %>% filter(pvalue < 0.05),
  "GSEA_CompuGenSets_MSigdb"=gseaC4_CellType1_kd1@result %>% filter(pvalue < 0.05),
  "GSEA_CellType_MSigdb"=gseaC8_CellType1_kd1@result %>% filter(pvalue < 0.05),
  "GSEA_BioProcess_GO"=gseBP_kd1@result %>% filter(pvalue < 0.05),
  "GSEA_MolFunc_GO"=gseMF_kd1@result %>% filter(pvalue < 0.05),
  "GSEA_CellComp_GO"=gseCC_kd1@result %>% filter(pvalue < 0.05)
)

openxlsx::write.xlsx(enrichTables_kd1,file="/tables/enrichmentAnalysisTables_CellType1_kd1sh.xlsx")

#gseaplot2(gseBP_CellType1kd2, geneSetID = 1:5)


#MsigDB in Enrichr:
#MSigDB_Computational, MSigDB_Oncogenic_Signatures,MSigDB_Hallmark_2020


# GSEA Cell Type 1 -----------------------------------------------------------
# Run this in Uppmax:
library(org.Hs.eg.db)
library(clusterProfiler)


sigGeneskd1df<- lfcShrink_mettl3kdkd1_ctrl %>% filter(padj <=0.05 & !is.na(padj)) %>% dplyr::select(gene,log2FoldChange) %>% distinct()
sigGeneslfckd1 <- sigGeneskd1df$log2FoldChange
names(sigGeneslfckd1) <- sigGeneskd1df$gene
sigGeneslfckd1 <- sort(sigGeneslfckd1,decreasing = T)

sigGeneskd2df<- lfcShrink_mettl3kdkd2_ctrl %>% filter(padj <=0.05 & !is.na(padj)) %>% dplyr::select(gene,log2FoldChange) %>% distinct()
sigGeneslfckd2 <- sigGeneskd2df$log2FoldChange
names(sigGeneslfckd2) <- sigGeneskd2df$gene
sigGeneslfckd2 <- sort(sigGeneslfckd2,decreasing = T)

fwrite(as.data.frame(sigGeneslfckd1),"/tables/sigGenesList_withLFC_kd1.tsv",row.names = T,sep="\t")

fwrite(as.data.frame(sigGeneslfckd2),"/tables/sigGenesList_withLFC_kd2.tsv",row.names = T,sep="\t")

# Run this in Uppmax:

gseBP_CellType1kd1<-gseGO(geneList = sigGeneslfckd1,OrgDb =org.Hs.eg.db,keyType = "SYMBOL",ont="BP",pvalueCutoff = 0.05)
gseBP_CellType1kd2<-gseGO(geneList = sigGeneslfckd2,OrgDb =org.Hs.eg.db,keyType = "SYMBOL",ont="BP",pvalueCutoff = 0.05)

load("/tables/gseaBP_CellType1_Prot1KD_kd1sh.Rda")
load("/tables/gseaBP_CellType1_Prot1KD_kd2sh.Rda")


# GSEA plots --------------------------------------------------------------

pdf("/plots/GSEA_cellType1_Control-sh_vs_Prot1KD_kd1_and_kd2.pdf", width = 15, height = 15)

gseaplot2(gseBP_CellType1kd1, geneSetID = 1:10,pvalue_table = T,title =str_wrap("GSEA top 10 terms - Cell Type 1: Significant genes Control-sh vs Prot1KD kd1",50) )

gseaplot2(gseBP_CellType1kd2, geneSetID = 1:10, title =str_wrap("GSEA top 10 terms - Cell Type 1 Significant genes Control-sh vs Prot1KD kd2",50))

heatplot(gseBP_CellType1kd1,showCategory = 70,foldChange = sigGeneslfckd1) + coord_flip() + labs(title="Heatmap GSEA: Genes and Terms (all) - GSEA Prot1KD kd1") + theme(axis.text.y = element_text(size=6)) + scale_fill_gradient2(low = "blue",mid = viridis(1),high = "salmon")

heatplot(gseBP_CellType1kd2,showCategory = 70,foldChange = sort(sigGeneslfckd2,decreasing = T)) + coord_flip() + labs(title="Heatmap GSEA: Genes and Terms (all) - GSEA Prot1KD kd2")  + theme(axis.text.y = element_text(size=6)) + scale_fill_gradient2(low = "blue",mid = viridis(1),high = "salmon")

dev.off()

gseaResults <- list(
  "GSEA_CellType1_kd1sh"=gseBP_CellType1kd1@result %>% filter(pvalue < 0.05) %>% arrange(desc(enrichmentScore),desc(rank)),
  "GSEA_CellType1_kd2sh"=gseBP_CellType1kd2@result %>% filter(pvalue < 0.05) %>% arrange(desc(enrichmentScore),desc(rank))
)

write.xlsx(gseaResults,file="/tables/GSEA_analysis_NCCday6_kd1sh_and_kd2sh.xlsx")

# Enrichment analysis NCC d6 ---------------------------------------------- 
#sigGeneslfckd1
goBP_kd1<-enrichGO(gene = names(sigGeneslfckd1),OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))

goBP_kd2<-enrichGO(gene = names(sigGeneslfckd2),OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))

pdf("/plots/enrichmentGOBiologicalProcess_kd1sh_kd2sh.pdf", width = 10, height = 14)

dotplot(goBP_kd1, orderBy="x",showCategory=30) + labs(title = "GO: Biological process \n Control-sh vs Prot1KD kd1-sh")
dotplot(goBP_kd2, orderBy="x", showCategory=30)  + labs(title = "GO: Biological process \n Control-sh vs Prot1KD kd2-sh")

dev.off()

goBPresults <- list(
  "goBP_kd1sh"=goBP_kd1,
  "goBP_kd2sh"=goBP_kd2
)

write.xlsx(goBPresults,file="/tables/enrichmentGOBiologicalProcess_CellType1_kd1sh_and_kd2sh.xlsx")




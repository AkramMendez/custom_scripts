library(DESeq2)
library(gprofiler2)
#Load counts data
args <- commandArgs(trailingOnly = TRUE)
featureCounts_file <- args[1]

featCounts <- fread(featureCounts_file, header = T,sep="\t", skip = 1)

colnames(featCounts) <- gsub("(\/.*\/|_S.*.ht2.srt.chrms.nodup.uniq.bam)","", colnames(featCounts),perl = T)

featCounts <- featCounts %>% mutate(Geneid=gsub("\\.\\d+","",featCounts$Geneid,perl = T))


genes <- gconvert(gsub("\\.\\d+","",featCounts$Geneid,perl = T),
         organism = "hsapiens",
         target = "HGNC",
         numeric_ns = "",
         mthreshold = Inf,
         filter_na = FALSE) %>% 
  dplyr::select(input,name)




featCounts <- left_join(featCounts,genes, by=c("Geneid"="input")) %>% 
  rename(gene=name) %>% 
  mutate(gene=if_else(is.na(gene),Geneid,gene)) %>% 
  filter(duplicated(gene)==F)

counts <- featCounts %>% dplyr::select(-c("gene","Geneid","Chr","Start","End","Strand","Length")) %>% 
  dplyr::select(contains("_Input")) %>% as.data.frame()


colnames(counts) <- gsub("Input_","R", colnames(counts))


rownames(counts) <- featCounts$gene

############ Calculate TPM from raw counts table for all genes and all samples ###################
# gtf_grch38<-"gencode.v36.chr_patch_hapl_scaff.basic.annotation.gtf"
# 
# txdb_grch38<-makeTxDbFromGFF(gtf_grch38,format = "gtf")
# 
# # Get transcript coordinates and select longest isoforms:
# txdf_hs<-AnnotationDbi::select(txdb_grch38,keys(txdb_grch38,"GENEID"),"TXNAME","GENEID")
# 
# 
# txdf_hs <- left_join(transcripts(txdb_grch38) %>% as.data.frame(), txdf_hs ,by=c("tx_name"="TXNAME")) %>% 
#   group_by(GENEID) %>% top_n(1,wt = width) %>% filter(., rank(desc(width), ties.method = "first")==1) %>% ungroup() %>% distinct()  
# 
# 
# txdf_hs$gene <- mapIds(org.Hs.eg.db,keys = gsub("\\..*","",txdf_hs$GENEID, perl = T) , column = "SYMBOL",keytype = "ENSEMBL",multiVals = "first")
# 
# txdf_hs <- txdf_hs %>% mutate(gene=if_else(is.na(gene),GENEID,gene))
# 
# geneLength_all_hg38.df <- txdf_hs %>% dplyr::select(gene,length=width)

geneLength_all_hg38.df <- featCounts %>% dplyr::select(gene,length=Length)

countssamples_TPM <- left_join(counts %>%
                               rownames_to_column("gene"),geneLength_all_hg38.df, by=c("gene"="gene")) %>% 
  pivot_longer(names_to = "sample",values_to = "counts",-c("gene","length")) %>%
  group_by(sample) %>% # Calculate values per each sample
  mutate(cpm= (counts/sum(as.numeric(counts))) * 10^6) %>% 
  mutate(rpk= counts/ (length/1000)) %>% # Calculate reads per kilobase dividing the counts by gene length, normalized to 1kb factor
  mutate(tpm=rpk/sum(rpk)*10^6) %>% ungroup()



##### LC day 4

metadata_cond1 <- NULL
metadata_cond1$sample<-as.factor(c("Cond1_Treat_R1","Cond1_Treat_R2","Cond1_Ctrl_R1","Cond1_Ctrl_R2"))
metadata_cond1$condition<-as.factor(c(rep("treatment",2),rep("control",2)))
metadata_cond1$condition<-relevel(metadata_cond1$condition,"control")
metadata_cond1<-as.data.frame(metadata_cond1)
rownames(metadata_cond1)<- c("Cond1_Treat_R1","Cond1_Treat_R2","Cond1_Ctrl_R1","Cond1_Ctrl_R2")


dds_Cond1 <- DESeqDataSetFromMatrix(countData = counts %>% dplyr::select(matches("Cond1")) , colData=metadata_cond1, design= ~ condition)
dds_Cond1<-DESeq(dds_Cond1,modelMatrixType = "standard", fitType = "local")
resultsNames(dds_Cond1)

normalized_counts_Cond1<- counts(dds_Cond1, normalized=TRUE)
normalized_counts_Cond1 <- normalized_counts_Cond1 %>% as.data.frame() %>%  mutate(gene=rownames(normalized_counts_Cond1))



res_Cond1 <- results(dds_Cond1, contrast = c("condition","treatment","control")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="Cond1Treat_vs_Cond1Control")

#Shrink LFC estimates:
lfcShrink_Cond1 <-lfcShrink(dds_Cond1,coef = "condition_treatment_vs_control",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="Cond1Treat_vs_Cond1Control")
lfcShrink_Cond1 <- lfcShrink_Cond1 %>% filter(!grepl("gene-",gene))


openxlsx::write.xlsx(list(
  "DiffExp_Cond1"=lfcShrink_Cond1 %>% dplyr::select(gene,baseMean,log2FoldChange,lfcSE,pvalue,padj,contrast) %>% filter(padj < 0.1 & !is.na(padj)) %>% arrange(desc(log2FoldChange),padj), 
  ), file = "diffExp_DEseq2__padj0.1_signif.xlsx",overwrite = TRUE)



library(EnhancedVolcano)

topsig_Cond1 <- lfcShrink_Cond1 %>% arrange(padj) %>% head(30) %>% dplyr::select(gene) %>% unlist()


EnhancedVolcano(lfcShrink_Cond1,lab = lfcShrink_Cond1$gene,x = "log2FoldChange",y = "padj",pCutoff = 0.05, FCcutoff = log2(1.2), pointSize = 0.1,labSize = 3.0, boxedLabels = F,maxoverlapsConnectors = 1000, drawConnectors = T,colConnectors = "lightgray",arrowheads = F,subtitle=str_wrap("LC day 4 \n Non-treatment vs Infected \n",50), title = "", legendPosition = 'right', legendLabSize = 11,legendIconSize = 3, axisLabSize = 11, caption = "", selectLab =c(topsig_Cond1,cov2genes,ciliaGenes))


######
#library(ggtext)
jointdotplotEnrichedTerms <- function(df1,df2,names,title){
  
  #file<-tolower(gsub(" ","",file))
  
  df1@result$GeneRatio<-sapply(df1@result$GeneRatio,function(x){eval(parse(text=x))})
  df2@result$GeneRatio<-sapply(df2@result$GeneRatio,function(x){eval(parse(text=x))})
  
  df1_filt<-df1@result %>% filter(pvalue < 0.05)  %>% mutate(condition=names[1], compose_score=-log10(pvalue)*GeneRatio) %>%  arrange(desc(compose_score)) %>% head(30)
  df2_filt<-df2@result %>% filter(pvalue < 0.05)  %>% mutate(condition=names[2],compose_score=-log10(pvalue)*GeneRatio) %>%  arrange(desc(compose_score)) %>% head(30)
  
  topTerms<-unique(union(df1_filt$Description,df2_filt$Description))
  commonTerms<-unique(intersect(df1_filt$Description,df2_filt$Description))
  diffTerms<-setdiff(topTerms,commonTerms)

  df<-rbind(df1_filt,df2_filt)
  df_filt<- df %>% filter(Description %in% topTerms) %>% arrange(desc(-log10(pvalue)),desc(GeneRatio))
  df_wide <-df_filt %>% pivot_wider(names_from = "condition",values_from = colnames(df_filt)[-c(1,2,10)])
  df_filt<- df_filt %>% arrange(dplyr::desc(compose_score))
  
  df_filt %>% ggplot(.,aes(condition,reorder(Description,compose_score), size=GeneRatio,color=-log10(pvalue))) + 
    geom_point(alpha=0.8) +
    scale_y_discrete(position = "right") +
    theme_light() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_markdown(angle=45,hjust = 1, vjust = 1,size=8),
          axis.text.y = element_markdown(size=6),aspect.ratio = 7) +
    scale_color_gradient(high="salmon",low = "blue") +
    labs(title=title) + xlab("") + ylab("") ->p
  return(p)
}



jointdotplotEnrichr <- function(df1,df2,df3,names,title){
  
  df1_filt<-df1 %>% filter(P.value < 0.05 )  %>% mutate(condition=names[1], compose_score= Combined.Score) %>%  arrange(desc(compose_score)) %>% head(30)
  df2_filt<-df2 %>% filter(P.value < 0.05 )  %>% mutate(condition=names[2],compose_score= Combined.Score) %>%  arrange(desc(compose_score)) %>% head(30)
  df3_filt<-df3 %>% filter(P.value < 0.05 ) %>% mutate(condition=names[3], compose_score= Combined.Score) %>%  arrange(desc(compose_score)) %>% head(30)
  %>% head(60)
  
  topTerms<-unique(union(union(df1_filt$Term,df2_filt$Term),df3_filt$Term))
  commonTerms<-unique(intersect(intersect(df1_filt$Term,df2_filt$Term),df3_filt$Term))
  diffTerms<-setdiff(topTerms,commonTerms)
  
  df<-rbind(df1_filt,df2_filt,df3_filt)
  df_filt<- df %>% filter(Term %in% topTerms) %>% arrange(desc(-log10(P.value)),desc(Odds.Ratio))
  df_wide <-df_filt %>% pivot_wider(names_from = "condition",values_from = colnames(df_filt)[-c(1,2,10)])
  df_filt<- df_filt %>% arrange(dplyr::desc(compose_score))
  
  
  df_filt %>% ggplot(.,aes(condition,reorder(Term,compose_score), size=Odds.Ratio,color=-log10(P.value))) + 
    geom_point(alpha=0.8) +
    scale_y_discrete(position = "right") +
    theme_light() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_markdown(angle=0,hjust = 0.5,size=8),
          axis.text.y = element_markdown(size=6),aspect.ratio = 5) +
    scale_color_gradient(high="salmon",low = "blue") +
    labs(title=title) + xlab("") + ylab("") ->p
  return(p)
}

############## Pathway analysis
library(clusterProfiler)
library(enrichR)

#filter(padj < 0.01 & abs(log2FoldChange) > 1) 
signif_Cond1 <- lfcShrink_Cond1 %>% filter(padj < 0.1 & !is.na(padj)

degs_Cond1 <- lfcShrink_Cond1 %>% filter(padj < 0.01 & !is.na(padj) &!(gene %in% cov2genes ) & abs(log2FoldChange) > 0.5849625) #0.5849625
#########################

signif_Cond1_entrez <- mapIds(org.Hs.eg.db,keys = signif_Cond1$gene, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")

goBP_Cond1<-enrichGO(gene = signif_Cond1$gene,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))

goMF_Cond1<-enrichGO(gene = signif_Cond1$gene,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))

goCC_Cond1<-enrichGO(gene = signif_Cond1$gene,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))

kegg_Cond1<-enrichKEGG(signif_Cond1_entrez,pvalueCutoff = 0.05)

kegg_Cond1 <- setReadable(kegg_Cond1, OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

reactome_Cond1<-enrichPathway(signif_Cond1_entrez,pvalueCutoff=0.05)
enrichrWP_Cond1<-enrichr(signif_Cond1$gene,"WikiPathways_2019_Human")

pdf("../figures/enrichmentAnalysis_significantGenes.pdf", width = 10, height = 8,onefile = T)

plotEnrich(enrichrWP_Cond1$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7)) + labs(title=str_wrap("Enrichr: WikiPathways \n Cond1 "),20)

dev.off()

goTables_Signif <- list(
  "GOBP_Cond1"=goBP_Cond1@result %>% filter(pvalue < 0.05),
  "GOMF_Cond1"=goMF_Cond1@result %>% filter(pvalue < 0.05),
  "GOCC_Cond1"=goCC_Cond1@result %>% filter(pvalue < 0.05),
  "WiPaths_Cond1"=enrichrWP_Cond1$WikiPathways_2019_Human %>% filter(P.value < 0.05),
  "KEGG_Cond1"=kegg_Cond1@result %>% filter(pvalue < 0.05),
  )

openxlsx::write.xlsx(goTables_Signif,"../tables/enrichmentAnalysisTables_significantGenes.xlsx", overwrite = T)
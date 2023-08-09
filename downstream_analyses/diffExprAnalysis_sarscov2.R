#!/usr/bin/env Rscript

library("Rsubread")
library("DESeq2")
library("tidyverse")
library("pheatmap")
library("RColorBrewer")
library("data.table")
library("ggtext")

args = commandArgs(trailingOnly=TRUE)

#inputdir <- "/crex/proj/nb_storage/mappings/covid19_infected_cells/nodup_uniq_alns/all"
#annotation_file<-"/crex/proj/nb_project/genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_chroms_concat.gtf"
inputdir <- "mappings/covid19_infected_cells/nodup_uniq_alns/all/"
annotation_file<-"/genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_chroms_concat.gtf"
#metadata_path <- args[2]
threads <- 1



setwd(inputdir)

bamfiles<-list.files(inputdir,pattern = ".chroms.nodup.uniq.sorted.bam$")
bamfiles<-grep("(Input)",bamfiles,perl = T,ignore.case = T,value = T)

#bamfiles2<-grep("(Input_2)",bamfiles,perl = T,ignore.case = T,value = T)

metadata<- NULL
metadata$sample<-as.factor(basename(bamfiles))
metadata$condition<-as.factor(c(rep("infected",6),rep("noninfected",2)))
metadata$condition<-relevel(metadata$condition,"noninfected")
#metadata$type<-as.factor(c("B.1","B.1.351","B.1.1.7","Vero"))
metadata$type<-as.factor(c("B.1","B.1","B.1.351","B.1.351","B.1.1.7","B.1.1.7","Vero","Vero"))
metadata$type<-relevel(metadata$type,"Vero")
metadata$run<-as.factor(rep(c("r1","r2"),4))
#metadata$group<-as.factor(paste(metadata$type,metadata$run,sep="_"))
#metadata$condition<-as.factor(c(rep("infected",3),rep("noninfected",1)))
#metadata$type<-as.factor(c("B.1","B.1.351","B.1.1.7","Vero"))
#metadata$run<-as.factor(rep(c("r1"),4))
metadata<-as.data.frame(metadata)
rownames(metadata)<-gsub("_S\\d+.*$","",basename(bamfiles),perl=TRUE)



#metadata<-fread(metadata_path, colClasses = "factor")

featureCounts_scov2<-featureCounts(bamfiles, annot.ext=annotation_file, isGTFAnnotationFile = TRUE, nthreads = threads)

counts_scov2<-featureCounts_scov2$counts
colnames(counts_scov2)<-gsub("_S\\d+.*$","",colnames(counts_scov2),perl=TRUE)
fwrite(as.data.frame(counts_scov2),"counts_infectedCellsAll_cov2_nodup_uniq.tsv",col.names=T,sep="\t", row.names =T)

#counts_scov2<-fread("mappings/covid19_infected_cells/nodup_uniq_alns/all/counts_infectedCellsAll_cov2_nodup_uniq.tsv",header = T,sep="\t")
#counts_scov2 <- counts_scov2 %>% column_to_rownames(.,"V1")
# Check that sample names match those from the design matrix metadata
if(!all(colnames(counts_scov2) %in% rownames(metadata))){
  stop("colnames(counts_scov2) not matching rownames(coldata_info)")
}




dds<-DESeqDataSetFromMatrix(countData = counts_scov2, colData=metadata, design= ~ type)

#design(dds)<- ~ type

dds<-DESeq(dds,modelMatrixType = "standard")

resultsNames(dds)

#results_all<-results(dds)
#wu_vs_vero<-results(dds,contrast = c("type","B.1","Vero"))

#_________________________________#
normalized_counts<- counts(dds, normalized=TRUE)

normalized_counts<-normalized_counts %>% as.data.frame() %>%  mutate(gene=rownames(normalized_counts))

fwrite(normalized_counts,"/plots/tables/normalizedCounts_diffexpAnalysis.tsv",sep="\t",col.names = T)
#results_all<-results(dds)

wu_vs_vero<-results(dds,contrast = c("type","B.1","Vero"),alpha = 0.01) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="wu_vs_vero")
uk_vs_vero<-results(dds,contrast = c("type","B.1.1.7","Vero")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="uk_vs_vero")
sa_vs_vero<-results(dds,contrast = c("type","B.1.351","Vero")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="sa_vs_vero")
# wu_vs_uk<-results(dds,contrast = c("type","B.1","B.1.1.7")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="wu_vs_uk")
# wu_vs_sa<-results(dds,contrast = c("type","B.1","B.1.351")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="wu_vs_sa")
# uk_vs_sa<-results(dds,contrast = c("type","B.1.1.7","B.1.351")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="uk_vs_sa")

#Shrink LFC estimates:
wu_vs_vero_lfcShrink<-lfcShrink(dds,coef = "type_B.1_vs_Vero",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="wu_vs_vero")
uk_vs_vero_lfcShrink<-lfcShrink(dds,coef = "type_B.1.1.7_vs_Vero",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="uk_vs_vero")
sa_vs_vero_lfcShrink<-lfcShrink(dds,coef = "type_B.1.351_vs_Vero",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="sa_vs_vero")
#wu_vs_uk_lfcShrhink<-lfcShrink(dds,coef = "type_WU_vs_UK",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="wu_vs_uk")
#wu_vs_sa_lfcShrhink<-lfcShrink(dds,coef = "type_WU_vs_SA",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="wu_vs_sa")
#uk_vs_sa_lfcShrhink<-lfcShrink(dds,coef = "type_UK_vs_SA",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="uk_vs_sa")



res_all<-rbind(wu_vs_vero,uk_vs_vero,sa_vs_vero)

res_all_lfcShrink<-rbind(wu_vs_vero_lfcShrink,uk_vs_vero_lfcShrink,sa_vs_vero_lfcShrink)

fwrite(res_all,"differentialExpression/diffexp_infected_vs_noninfected_by_strain.tsv",col.names = T,sep="\t")
fwrite(res_all_lfcShrink,"differentialExpression/diffexp_infected_vs_noninfected_by_strain_lfcShrink.tsv",col.names = T,sep="\t")
#res_all_lfcShrink<-fread("differentialExpression/diffexp_infected_vs_noninfected_by_strain_lfcShrink.tsv",header = T,sep="\t")

res_all_lfcShrink %>% filter(!(gene %in% cov2genes) & !grepl("gene-",gene)) %>% arrange(dplyr::desc(contrast),dplyr::desc(log2FoldChange),padj) %>% dplyr::select(contrast,gene,symbol,gene_name,log2FoldChange,pvalue,padj,lfcSE,baseMean) %>% distinct() %>% 
  openxlsx::write.xlsx(.,"/plots/tables/Table_differentialExpressionAnalysis.xls")

deg_tables <- list(
  "B.1_vs_NonInfectedVero"= res_all_lfcShrink %>% filter(!(gene %in% cov2genes) & !grepl("gene-",gene) & contrast=="wu_vs_vero" & abs(log2FoldChange) >1 & padj < 0.01 &! is.na(padj)) %>% arrange(dplyr::desc(log2FoldChange),padj) %>%
  dplyr::select(gene,symbol,gene_name,log2FoldChange,pvalue,padj,lfcSE,baseMean) %>% distinct() ,
  "B.1.1.7_vs_NonInfectedVero"= res_all_lfcShrink %>% filter(!(gene %in% cov2genes) & !grepl("gene-",gene) & contrast=="uk_vs_vero" & abs(log2FoldChange) >1 & padj < 0.01 &! is.na(padj)) %>% arrange(dplyr::desc(log2FoldChange),padj) %>%
    dplyr::select(gene,symbol,gene_name,log2FoldChange,pvalue,padj,lfcSE,baseMean) %>% distinct(),
  "B.1.351_vs_NonInfectedVero"= res_all_lfcShrink %>% filter(!(gene %in% cov2genes) & !grepl("gene-",gene) & contrast=="sa_vs_vero" & abs(log2FoldChange) >1 & padj < 0.01 &! is.na(padj)) %>% arrange(dplyr::desc(log2FoldChange),padj) %>%
    dplyr::select(gene,symbol,gene_name,log2FoldChange,pvalue,padj,lfcSE,baseMean) %>% distinct()
)

openxlsx::write.xlsx(deg_tables,"/plots/tables/DiffExprGenes_by_strain.xls")

cov2genes<-c("S","ORF8","ORF7a","ORF7b","ORF6","ORF3a","ORF1ab","ORF10","N","M","E","ORF1a")

res_all_lfcShrink<-res_all_lfcShrink %>% filter(!(gene %in% cov2genes))

library("AnnotationDbi")
library("org.Hs.eg.db")
library(pathview)
library(gage)
library(gageData)
library(GO.db)

# set up kegg database
kg.hsa <- kegg.gsets(species="hsa")
kegg.sigmet.gs <- kg.hsa$kg.sets[kg.hsa$sigmet.idx]
kegg.dise.gs <- kg.hsa$kg.sets[kg.hsa$dise.idx]

# set up go database
go.hs <- go.gsets(species="human")
go.bp.gs <- go.hs$go.sets[go.hs$go.subs$BP]
go.mf.gs <- go.hs$go.sets[go.hs$go.subs$MF]
go.cc.gs <- go.hs$go.sets[go.hs$go.subs$CC]

res_all_lfcShrink$symbol <- mapIds(org.Hs.eg.db,keys = res_all_lfcShrink$gene, column = "SYMBOL",keytype = "SYMBOL",multiVals = "first")

res_all_lfcShrink$entrez <- mapIds(org.Hs.eg.db,keys = res_all_lfcShrink$gene, column = "ENTREZID",keytype = "SYMBOL",multiVals="list")
res_all_lfcShrink$gene_name <- mapIds(org.Hs.eg.db,keys = res_all_lfcShrink$gene, column = "GENENAME",keytype = "SYMBOL",multiVals = "first")


 
logfc_wu_vs_vero<- res_all_lfcShrink %>% filter(contrast=="wu_vs_vero" &! grepl("gene-",gene))
logfc_uk_vs_vero<- res_all_lfcShrink %>% filter(contrast=="uk_vs_vero" &! grepl("gene-",gene))
logfc_sa_vs_vero<- res_all_lfcShrink %>% filter(contrast=="sa_vs_vero" &! grepl("gene-",gene))

foldchanges_wu_vs_vero<-logfc_wu_vs_vero$log2FoldChange
foldchanges_uk_vs_vero<-logfc_uk_vs_vero$log2FoldChange
foldchanges_sa_vs_vero<-logfc_sa_vs_vero$log2FoldChange

names(foldchanges_wu_vs_vero)<-logfc_wu_vs_vero$entrez
names(foldchanges_uk_vs_vero)<-logfc_uk_vs_vero$entrez
names(foldchanges_sa_vs_vero)<-logfc_sa_vs_vero$entrez


##############################################
fc.kegg.sigmet_wu_vs_vero<- gage(foldchanges_wu_vs_vero,gsets = kegg.sigmet.gs)
fc.kegg.sigmet_uk_vs_vero<- gage(foldchanges_uk_vs_vero,gsets = kegg.sigmet.gs)
fc.kegg.sigmet_sa_vs_vero<- gage(foldchanges_sa_vs_vero,gsets = kegg.sigmet.gs)


fc.kegg.disease_wu_vs_vero<-gage(foldchanges_wu_vs_vero,gsets = kegg.dise.gs)
fc.kegg.disease_uk_vs_vero<-gage(foldchanges_uk_vs_vero,gsets = kegg.dise.gs)
fc.kegg.disease_sa_vs_vero<-gage(foldchanges_sa_vs_vero,gsets = kegg.dise.gs)

fc.go.bp_wu_vs_vero<-gage(foldchanges_wu_vs_vero,gsets = go.bp.gs)
fc.go.bp_uk_vs_vero<-gage(foldchanges_uk_vs_vero,gsets = go.bp.gs)
fc.go.bp_sa_vs_vero<-gage(foldchanges_sa_vs_vero,gsets = go.bp.gs)

fc.go.mf_wu_vs_vero<-gage(foldchanges_wu_vs_vero,gsets = go.mf.gs)
fc.go.mf_uk_vs_vero<-gage(foldchanges_uk_vs_vero,gsets = go.mf.gs)
fc.go.mf_sa_vs_vero<-gage(foldchanges_sa_vs_vero,gsets = go.mf.gs)

fc.cc.mf_wu_vs_vero<-gage(foldchanges_wu_vs_vero,gsets = go.cc.gs)
fc.cc.mf_uk_vs_vero<-gage(foldchanges_uk_vs_vero,gsets = go.cc.gs)
fc.cc.mf_sa_vs_vero<-gage(foldchanges_sa_vs_vero,gsets = go.cc.gs)

#####################################################
#convert to data frame:
# fc.signal.less_wu_vs_vero<-as.data.frame(fc.kegg.sigmet_wu_vs_vero$less) %>% filter(greater.p.val<0.05) %>% mutate(contrast="wu_vs_vero",enrichment="less")
# fc.signal.less_uk_vs_vero<-as.data.frame(fc.kegg.sigmet_uk_vs_vero$less) %>% filter(p.val<0.05) %>% mutate(contrast="uk_vs_vero",enrichment="less")
# fc.signal.less_sa_vs_vero<-as.data.frame(fc.kegg.sigmet_sa_vs_vero$less) %>% filter(p.val<0.05) %>% mutate(contrast="sa_vs_vero",enrichment="less")
# 
# fc.signal.greater_wu_vs_vero<-as.data.frame(fc.kegg.sigmet_wu_vs_vero$greater) %>% filter(p.val<0.05) %>% mutate(contrast="wu_vs_vero",enrichment="greater")
# fc.signal.greater_uk_vs_vero<-as.data.frame(fc.kegg.sigmet_uk_vs_vero$greater) %>% filter(p.val<0.05) %>% mutate(contrast="uk_vs_vero",enrichment="greater")
# fc.signal.greater_sa_vs_vero<-as.data.frame(fc.kegg.sigmet_sa_vs_vero$greater) %>% filter(p.val<0.05) %>% mutate(contrast="sa_vs_vero",enrichment="greater")
# 
# fwrite(rbind(fc.signal.less_wu_vs_vero,fc.signal.less_uk_vs_vero),"differentialExpression/pathway_analysis/signal_enrichment_foldChanges_less.tsv", col.names=T,sep="\t")
# 
# fwrite(rbind(fc.signal.greater_wu_vs_vero,fc.signal.greater_uk_vs_vero,fc.signal.greater_sa_vs_vero),"differentialExpression/pathway_analysis/signal_enrichment_foldChanges_greater.tsv", col.names=T,sep="\t")

setwd("/differentialExpression/pathway_analysis/plots/")
fc.kegg.disease_wu_vs_vero$greater %>% as.data.frame() %>% filter(q.val < 0.05 &! is.na(q.val)) %>% mutate(name=rownames(.)) %>% head(50) %>% ggplot(.,aes(reorder(name,dplyr::desc(q.val)),set.size)) + geom_bar(stat = "identity", fill="steelblue") + coord_flip() + theme_light() + ylab("Set size") + xlab("Kegg Disease") + labs(title = "non-infected vs infected A.1",subtitle="Enrichment analysis (Kegg disease)") -> p1
ggsave(filename = "kegg_disease_vero_vs_wu_fdr0.05.pdf",plot = p1,height = 8,width = 8,dpi=300,device = "pdf")

fc.kegg.sigmet_wu_vs_vero$greater %>% as.data.frame() %>% filter(q.val < 0.05 &! is.na(q.val)) %>% mutate(name=rownames(.)) %>% head(50) %>% ggplot(.,aes(reorder(name,dplyr::desc(q.val)),set.size)) + geom_bar(stat = "identity", fill="steelblue") + coord_flip() + theme_light() + ylab("Set size") + xlab("Kegg Signaling") + labs(title = "non-infected vs infected A.1",subtitle="Enrichment analysis (Kegg signaling)")  -> p2
ggsave(filename = "kegg_signaling_vero_vs_wu_fdr0.05.pdf",plot = p2,height = 8,width = 8,dpi=300,device = "pdf")

fc.go.mf_wu_vs_vero$greater %>% as.data.frame() %>% filter(q.val < 0.05 &! is.na(q.val)) %>% mutate(name=rownames(.)) %>% head(50) %>% ggplot(.,aes(reorder(name,dplyr::desc(q.val)),set.size)) + geom_bar(stat = "identity", fill="steelblue") + coord_flip() + theme_light() + ylab("Set size") + xlab("Molecular function") + labs(title = "non-infected vs infected A.1",subtitle="Enrichment analysis (Molecular function)") -> p3
ggsave(filename = "molecular_function_vero_vs_wu_fdr0.05.pdf",plot = p3,height = 8,width = 8,dpi=300,device = "pdf")

fc.go.bp_wu_vs_vero$greater %>% as.data.frame() %>% filter(q.val < 0.05 &! is.na(q.val)) %>% mutate(name=rownames(.)) %>% head(50) %>% ggplot(.,aes(reorder(name,dplyr::desc(q.val)),set.size)) + geom_bar(stat = "identity", fill="steelblue") + coord_flip() + theme_light() + ylab("Set size") + xlab("Biological process") + labs(title = "non-infected vs infected A.1",subtitle="Enrichment analysis (Biological process)") -> p4
ggsave(filename = "biological_process_vero_vs_wu_fdr0.05.pdf",plot = p4,height = 8,width = 8,dpi=300,device = "pdf")

fc.kegg.disease_uk_vs_vero$greater %>% as.data.frame() %>% mutate(name=rownames(.)) %>% head(50) %>% ggplot(.,aes(reorder(name,dplyr::desc(q.val)),set.size)) + geom_bar(stat = "identity", fill="steelblue") + coord_flip() + theme_light() + ylab("Set size") + xlab("Kegg Disease") + labs(title = "non-infected vs infected B.1.1.7",subtitle="Enrichment analysis (Kegg disease)") -> p5
ggsave(filename = "kegg_disease_vero_vs_uk_fdr0.05.pdf",plot = p5,height = 8,width = 8,dpi=300,device = "pdf")

fc.kegg.sigmet_uk_vs_vero$greater %>% as.data.frame() %>% mutate(name=rownames(.)) %>% head(50) %>% ggplot(.,aes(reorder(name,dplyr::desc(q.val)),set.size)) + geom_bar(stat = "identity", fill="steelblue") + coord_flip() + theme_light() + ylab("Set size") + xlab("Kegg Signaling") + labs(title = "non-infected vs infected B.1.1.7",subtitle="Enrichment analysis (Kegg signaling)") -> p6
ggsave(filename = "kegg_signaling_vero_vs_uk_fdr0.05.pdf",plot = p6,height = 8,width = 8,dpi=300,device = "pdf")

fc.go.mf_uk_vs_vero$greater %>% as.data.frame() %>% mutate(name=rownames(.)) %>% head(50) %>% ggplot(.,aes(reorder(name,dplyr::desc(q.val)),set.size)) + geom_bar(stat = "identity", fill="steelblue") + coord_flip() + theme_light() + ylab("Set size") + xlab("Molecular function") + labs(title = "non-infected vs infected B.1.1.7",subtitle="Enrichment analysis (Molecular function)")-> p7
ggsave(filename = "molecular_function_vero_vs_uk_fdr0.05.pdf",plot = p7,height = 8,width = 8,dpi=300,device = "pdf")

fc.go.bp_uk_vs_vero$greater %>% as.data.frame() %>% mutate(name=rownames(.)) %>% head(50) %>% ggplot(.,aes(reorder(name,dplyr::desc(q.val)),set.size)) + geom_bar(stat = "identity", fill="steelblue") + coord_flip() + theme_light() + ylab("Set size") + xlab("Biological process") + labs(title = "non-infected vs infected B.1.1.7",subtitle="Enrichment analysis (Biological process)")-> p8
ggsave(filename = "biological_process_vero_vs_uk_fdr0.05.pdf",plot = p8,height = 8,width = 8,dpi=300,device = "pdf")

fc.kegg.disease_sa_vs_vero$greater %>% as.data.frame() %>% mutate(name=rownames(.)) %>% head(50) %>% ggplot(.,aes(reorder(name,dplyr::desc(q.val)),set.size)) + geom_bar(stat = "identity", fill="steelblue") + coord_flip() + theme_light() + ylab("Set size") + xlab("Kegg Disease") + labs(title = "non-infected vs infected B.1.351",subtitle="Enrichment analysis (Kegg disease)") -> p9
ggsave(filename = "kegg_disease_vero_vs_sa_fdr0.05.pdf",plot = p9,height = 8,width = 8,dpi=300,device = "pdf")

fc.kegg.sigmet_sa_vs_vero$greater %>% as.data.frame() %>% mutate(name=rownames(.)) %>% head(50) %>% ggplot(.,aes(reorder(name,dplyr::desc(q.val)),set.size)) + geom_bar(stat = "identity", fill="steelblue") + coord_flip() + theme_light() + ylab("Set size") + xlab("Kegg Signaling") + labs(title = "non-infected vs infected B.1.351",subtitle="Enrichment analysis (Kegg signaling)")-> p10
ggsave(filename = "kegg_signaling_vero_vs_sa_fdr0.05.pdf",plot = p10,height = 8,width = 8,dpi=300,device = "pdf")

fc.go.mf_sa_vs_vero$greater %>% as.data.frame() %>% mutate(name=rownames(.)) %>% head(50) %>% ggplot(.,aes(reorder(name,dplyr::desc(q.val)),set.size)) + geom_bar(stat = "identity", fill="steelblue") + coord_flip() + theme_light() + ylab("Set size") + xlab("Molecular function") + labs(title = "non-infected vs infected B.1.351",subtitle="Enrichment analysis (Molecular function)")->p11
ggsave(filename = "molecular_function_vero_vs_sa_fdr0.05.pdf",plot = p11,height = 8,width = 8,dpi=300,device = "pdf")

fc.go.bp_sa_vs_vero$greater %>% as.data.frame() %>% mutate(name=rownames(.)) %>% head(50) %>% ggplot(.,aes(reorder(name,dplyr::desc(q.val)),set.size)) + geom_bar(stat = "identity", fill="steelblue") + coord_flip() + theme_light() + ylab("Set size") + xlab("Biological process") + labs(title = "non-infected vs infected B.1.351",subtitle="Enrichment analysis (Biological process)")-> p12
ggsave(filename = "biological_process_vero_vs_sa_fdr0.05.pdf",plot = p12,height = 8,width = 8,dpi=300,device = "pdf")
#####################################
#Dotplot 
#Enrichment analysis with clusterProfiler:
library(clusterProfiler)
sig_wu<-logfc_wu_vs_vero %>% filter(padj<0.01 &! is.na(padj) & abs(log2FoldChange)>1 &!(gene %in% cov2genes ))
sig_wu<-sig_wu$symbol
names(sig_wu)<-sig_wu
sig_wu_entrez<-mapIds(org.Hs.eg.db,keys = sig_wu, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
names(sig_wu_entrez)<-sig_wu_entrez

sig_uk<-logfc_uk_vs_vero %>% filter(padj<0.01 &! is.na(padj) & abs(log2FoldChange)>1 &!(gene %in% cov2genes ))
sig_uk<-sig_uk$symbol
names(sig_uk)<-sig_wu
#background_uk<-unique(logfc_uk_vs_vero$symbol)
#names(background_uk)<-background_uk
sig_uk_entrez<-mapIds(org.Hs.eg.db,keys = sig_uk, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
names(sig_uk_entrez)<-sig_uk_entrez
#go_bp_uk %>% slot("result") %>% as.tibble()

sig_sa<-logfc_sa_vs_vero %>% filter(padj<0.01 &! is.na(padj) & abs(log2FoldChange)>1 &!(gene %in% cov2genes ))
sig_sa<-sig_sa$symbol
names(sig_sa)<-sig_sa
background_sa<-unique(logfc_sa_vs_vero$symbol)
names(background_sa)<-background_sa
sig_sa_entrez<-mapIds(org.Hs.eg.db,keys = sig_sa, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
names(sig_sa_entrez)<-sig_sa_entrez

go_bp_wu<-enrichGO(gene = sig_wu,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
go_bp_uk<-enrichGO(gene = sig_uk,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
go_bp_sa<-enrichGO(gene = sig_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))

#Gene set enrichment analysis
sig_wu_gse<-logfc_wu_vs_vero %>% filter(padj<0.01 &! is.na(padj) & abs(log2FoldChange)>1 &!(gene %in% cov2genes )) %>% dplyr::select(log2FoldChange) %>% unlist()
names(sig_wu_gse)<-logfc_wu_vs_vero %>% filter(padj<0.01 &! is.na(padj) & abs(log2FoldChange)>1 &!(gene %in% cov2genes )) %>% dplyr::select(gene) %>% unlist()
sig_wu_gse<-sort(sig_wu_gse,decreasing = T)
gse_bp_wu<-gseGO(geneList = sig_wu_gse,OrgDb =org.Hs.eg.db,keyType = "SYMBOL",ont="BP",pvalueCutoff = 0.05)
cnetplot(gse_bp_wu,foldChange = sig_wu_gse,cex_gene=0.1,showCategory = 5,categorySize="pvalue",cex_label_gene=0.5,cex_label_category=0.5)

go_mf_wu<-enrichGO(gene = sig_wu,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
go_mf_uk<-enrichGO(gene = sig_uk,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
go_mf_sa<-enrichGO(gene = sig_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))

go_cc_wu<-enrichGO(gene = sig_wu,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
go_cc_uk<-enrichGO(gene = sig_uk,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
go_cc_sa<-enrichGO(gene = sig_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))

kegg_wu<-enrichKEGG(sig_wu_entrez,pvalueCutoff = 0.05)
kegg_uk<-enrichKEGG(sig_uk_entrez,pvalueCutoff = 0.05)
kegg_sa<-enrichKEGG(sig_sa_entrez,pvalueCutoff = 0.05)

kegg_wu <- setReadable(kegg_wu, OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
kegg_uk <- setReadable(kegg_uk, OrgDb = org.Hs.eg.db,keyType = "ENTREZID")
kegg_sa <- setReadable(kegg_sa, OrgDb = org.Hs.eg.db,keyType = "ENTREZID")

reactome_wu<-enrichPathway(sig_wu_entrez,pvalueCutoff=0.05)
reactome_uk<-enrichPathway(sig_uk_entrez,pvalueCutoff=0.05)
reactome_sa<-enrichPathway(sig_sa_entrez,pvalueCutoff=0.05)

enrichrWP_wu<-enrichr(sig_wu,"WikiPathways_2019_Human")
enrichrWP_uk<-enrichr(sig_uk,"WikiPathways_2019_Human")
enrichrWP_sa<-enrichr(sig_sa,"WikiPathways_2019_Human")

enrichrCov2sets_wu<-enrichr(sig_wu,"COVID-19_Related_Gene_Sets")
enrichrCov2sets_uk<-enrichr(sig_uk,"COVID-19_Related_Gene_Sets")
enrichrCov2sets_sa<-enrichr(sig_sa,"COVID-19_Related_Gene_Sets")


#go_bp_sa %>% slot("result") %>% as.tibble()


# pdf("/plots/biological_process_differentiallyExpressedGenes_all.pdf", width = 10, height = 8,onefile = T)
# mydotplot(go_bp_wu,showCategory=30,x="Count") + labs(title = "GO:Biological process \n Differentially expressed genes Vero vs WU")
# mydotplot(go_bp_uk,showCategory=30,x="Count") + labs(title = "GO:Biological process \n Differentially expressed genes Vero vs UK")
# mydotplot(go_bp_sa,showCategory=30,x="Count") + labs(title = "GO:Biological process \n Differentially expressed genes Vero vs SA")
# dev.off()
tables_outdir="/plots/tables/"

pdf("/plots/compiled_figures/enrichmentGO_topTerms_DEGs_by_Strain_jointPlot.pdf", width = 10, height = 8,onefile = T)

jointdotplot(go_bp_wu,go_bp_uk,go_bp_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "GO:Biological process \n DEGs B.1, B.1.1.7, B.1.351")
jointdotplot(go_mf_wu,go_mf_uk,go_mf_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "GO:Molecular function \n DEGs B.1, B.1.1.7, B.1.351")
jointdotplot(go_cc_wu,go_cc_uk,go_cc_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "GO:Cellular component \n DEGs B.1, B.1.1.7, B.1.351")
jointdotplot(kegg_wu,kegg_uk,kegg_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "Kegg \n All DEGs",outdir = tables_outdir, file="kegg_allDEGs" )

plotEnrich(enrichrCov2sets_wu$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + labs(title="Enrichr: COVID-19 Related Gene Sets \ B.1 ")
plotEnrich(enrichrCov2sets_uk$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + labs(title="Enrichr: COVID-19 Related Gene Sets \ B.1.1.7 ")
plotEnrich(enrichrCov2sets_sa$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + labs(title="Enrichr: COVID-19 Related Gene Sets \ B.1.351 ")

plotEnrich(enrichrWP_wu$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7)) + labs(title="Enrichr: WikiPathways \ B.1")
plotEnrich(enrichrWP_uk$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7)) + labs(title="Enrichr: WikiPathways \ B.1.1.7")
plotEnrich(enrichrWP_sa$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7)) + labs(title="Enrichr: WikiPathways \ B.1.351")

dev.off()

goTables_DEGs_all<-list("BP_DEGs_B.1"=go_bp_wu@result %>% filter(pvalue < 0.05) ,"BP_DEGs_B.1.1.7"=go_bp_uk@result %>% filter(pvalue < 0.05) ,"BP_DEGs_B.1.351"=go_bp_sa@result %>% filter(pvalue < 0.05) ,"MF_DEGs_B.1"=go_mf_wu@result %>% filter(pvalue < 0.05) ,"MF_DEGs_B.1.1.7"=go_mf_uk@result %>% filter(pvalue < 0.05) ,"MF_DEGs_B.1.351"=go_mf_sa@result %>% filter(pvalue < 0.05) ,"CC_DEGs_B.1"=go_cc_wu@result %>% filter(pvalue < 0.05) ,"CC_DEGs_B.1.1.7"=go_cc_uk@result %>% filter(pvalue < 0.05) ,"CC_DEGs_B.1.351"=go_cc_sa@result %>% filter(pvalue < 0.05) ,"KEGG_DEGs_B.1"=kegg_wu@result %>% filter(pvalue < 0.05) ,"KEGG_DEGs_B.1.1.7"=kegg_uk@result %>% filter(pvalue < 0.05),"KEGG_DEGs_B.1.351"=kegg_sa@result %>% filter(pvalue < 0.05),"EnrichrCov2Sets_B.1"=enrichrCov2sets_wu$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"EnrichrCov2Sets_B.1.1.7"=enrichrCov2sets_uk$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"EnrichrCov2Sets_B.1.351"=enrichrCov2sets_sa$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"WikiPaths_B.1"=enrichrWP_wu$WikiPathways_2019_Humangenes,"WikiPaths_DEGs_B.1.1.7"=enrichrWP_wu$WikiPathways_2019_Humangenes,"WikiPaths_DEGs_B.1.351"=enrichrWP_wu$WikiPathways_2019_Human %>% filter(P.value < 0.05))
openxlsx::write.xlsx(goTables_DEGs_all,"/plots/compiled_figures/tables/EnrcihmentAnalysis_DEGs_by_Strain.xlsx")
#goTables_DEGs_all<-openxlsx::read.xlsx("/plots/compiled_figures/tables/EnrcihmenAnalysis_DEGs_by_Strain.xlsx")
#common_diffexp_all<-res_all_lfcShrink

diffexpGenes_wu<-logfc_wu_vs_vero %>% filter(abs(log2FoldChange) > 1 &!(gene %in% cov2genes) & padj < 0.01 &! is.na(padj)) %>% dplyr::select(gene) %>% distinct() %>% unlist() 
diffexpGenes_uk<-logfc_uk_vs_vero %>% filter(abs(log2FoldChange) > 1 &!(gene %in% cov2genes)  & padj < 0.01 &! is.na(padj)) %>% dplyr::select(gene) %>% distinct() %>% unlist() 
diffexpGenes_sa<-logfc_sa_vs_vero %>% filter(abs(log2FoldChange) > 1 &!(gene %in% cov2genes)  & padj < 0.01 &! is.na(padj)) %>% dplyr::select(gene) %>% distinct() %>% unlist()

diffexpGenesUp_wu<-logfc_wu_vs_vero %>% filter(log2FoldChange > 1 &!(gene %in% cov2genes)  & padj < 0.01 &! is.na(padj)) %>% dplyr::select(gene) %>% distinct() %>% unlist() 
diffexpGenesUp_uk<-logfc_uk_vs_vero %>% filter(log2FoldChange > 1 &!(gene %in% cov2genes) & padj < 0.01 &! is.na(padj)) %>% dplyr::select(gene) %>% distinct() %>% unlist() 
diffexpGenesUp_sa<-logfc_sa_vs_vero %>% filter(log2FoldChange > 1 &!(gene %in% cov2genes) & padj < 0.01 &! is.na(padj)) %>% dplyr::select(gene) %>% distinct() %>% unlist() 

diffexpGenesUp_wu_entrez<-mapIds(org.Hs.eg.db,keys = diffexpGenesUp_wu, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
diffexpGenesUp_uk_entrez<-mapIds(org.Hs.eg.db,keys = diffexpGenesUp_uk, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
diffexpGenesUp_sa_entrez<-mapIds(org.Hs.eg.db,keys = diffexpGenesUp_sa, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")


diffexpGenesDown_wu<-logfc_wu_vs_vero %>% filter(log2FoldChange < -1 &!(gene %in% cov2genes)  & padj < 0.01 &! is.na(padj)) %>% dplyr::select(gene) %>% distinct() %>% unlist() 
diffexpGenesDown_uk<-logfc_uk_vs_vero %>% filter(log2FoldChange < -1 &!(gene %in% cov2genes)  & padj < 0.01 &! is.na(padj)) %>% dplyr::select(gene) %>% distinct() %>% unlist() 
diffexpGenesDown_sa<-logfc_sa_vs_vero %>% filter(log2FoldChange < -1 &!(gene %in% cov2genes)  & padj < 0.01 &! is.na(padj)) %>% dplyr::select(gene) %>% distinct() %>% unlist() 

diffexpGenesDown_wu_entrez<-mapIds(org.Hs.eg.db,keys = diffexpGenesDown_wu, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
diffexpGenesDown_uk_entrez<-mapIds(org.Hs.eg.db,keys = diffexpGenesDown_uk, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
diffexpGenesDown_sa_entrez<-mapIds(org.Hs.eg.db,keys = diffexpGenesDown_sa, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")

common_diffexp_genes<-intersect(intersect(diffexpGenes_wu,diffexpGenes_uk),diffexpGenes_sa)
fwrite(as.data.frame(common_diffexp_genes),"/plots/compiled_figures/tables/common_DEGs_all.tsv", col.names = F, sep="\t")
common_diffexp_genes_entrez<-mapIds(org.Hs.eg.db,keys = common_diffexp_genes, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")

goBP_commonDiffExpGenes<-enrichGO(gene=common_diffexp_genes,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F)%>% gofilter(.,level = c(7:10))
goMF_commonDiffExpGenes<-enrichGO(gene=common_diffexp_genes,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goCC_commonDiffExpGenes<-enrichGO(gene=common_diffexp_genes,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
kegg_commonDiffExpGenes<-enrichKEGG(common_diffexp_genes_entrez,pvalueCutoff = 0.05)
reactome_commonDiffExpGenes<-enrichPathway(common_diffexp_genes_entrez,pvalueCutoff=0.05)
enrichrWP_commonDiffExpGenes<-enrichr(common_diffexp_genes,"WikiPathways_2019_Human")
enrichrCov2sets_commonDiffExpGenes<-enrichr(common_diffexp_genes,"COVID-19_Related_Gene_Sets")

pdf("/plots/compiled_figures/enrichmentGO_topTerms_commonDEGs.pdf", width = 12, height = 16)
 mydotplot(goBP_commonDiffExpGenes,showCategory=30,x="Count", title="GO:Biological Process \n Common differentially expressed genes \n B.1, B.1.1.7, B.1.351")
 mydotplot(goMF_commonDiffExpGenes,showCategory=30,x="Count",title="GO:Molecular function \n Common differentially expressed genes \n B.1, B.1.1.7, B.1.351")
 mydotplot(goCC_commonDiffExpGenes,showCategory=30,x="Count", title="GO:Cellular component \n Common differentially expressed genes \n B.1, B.1.1.7, B.1.351")
 mydotplot(kegg_commonDiffExpGenes,showCategory=30,x="Count", title="KEGG pathways \n Common differentially expressed genes \n B.1, B.1.1.7, B.1.351")

# jointdotplot(goBP_commonDiffExpGenes,names = c("WU & UK & SA"),title = "GO:Biological process \n Common DEGs")
# jointdotplot(goMF_commonDiffExpGenes,names =  c("WU & UK & SA"),title = "GO:Molecular function \n Common DEGs")
# jointdotplot(goCC_commonDiffExpGenes,names =  c("WU & UK & SA"),title = "GO:Cellular component \n Common DEGs")
# jointdotplot(kegg_commonDiffExpGenes,names =  c("WU & UK & SA"),title = "Kegg \n Common DEGs")

plotEnrich(enrichrCov2sets_commonDiffExpGenes$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7)) + labs(title="Enrichr: COVID-19 Related Gene Sets \ Common differentially expressed genes ")
plotEnrich(enrichrWP_commonDiffExpGenes$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))   + labs(title="Enrichr: WikiPathways \ Common differentially expressed genes ")
dev.off()

goTables_common_DEGS<-list("BP_commonDEGs"=goBP_commonDiffExpGenes@result %>% filter(pvalue < 0.05),"MF_commonDEGs"=goMF_commonDiffExpGenes@result %>% filter(pvalue < 0.05),"CC_commonDEGs"=goCC_commonDiffExpGenes@result %>% filter(pvalue < 0.05),"KEGG_common_DEGs"=kegg_commonDiffExpGenes,"EnrichrCov2Sets_commonDEGs"=enrichrCov2sets_commonDiffExpGenes$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"WikiPath_commonDEgs"=enrichrWP_commonDiffExpGenes$WikiPathways_2019_Human %>% filter(P.value < 0.05))
openxlsx::write.xlsx(list("common_DEGs"=common_diffexp_genes),"/plots/compiled_figures/tables/commonDEGs_geneName_list.xlsx")
openxlsx::write.xlsx(goTables_common_DEGS,"/plots/compiled_figures/tables/EnrcihmenAnalysis_commonDEGs_All.xlsx")

#Enrichment Up regulated genes by condition
goBP_diffexpGenesUp_wu<-enrichGO(gene=diffexpGenesUp_wu,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goMF_diffexpGenesUp_wu<-enrichGO(gene=diffexpGenesUp_wu,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goCC_diffexpGenesUp_wu<-enrichGO(gene=diffexpGenesUp_wu,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
kegg_diffexpGenesUp_wu<-enrichKEGG(diffexpGenesUp_wu_entrez,pvalueCutoff = 0.05)
reactome_diffexpGenesUp_wu<-enrichPathway(diffexpGenesUp_wu_entrez,pvalueCutoff=0.05)
enrichrWP_diffexpGenesUp_wu<-enrichr(diffexpGenesUp_wu,"WikiPathways_2019_Human")
enrichrCov2sets_diffexpGenesUp_wu<-enrichr(diffexpGenesUp_wu,"COVID-19_Related_Gene_Sets")

goBP_diffexpGenesUp_uk<-enrichGO(gene=diffexpGenesUp_uk,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goMF_diffexpGenesUp_uk<-enrichGO(gene=diffexpGenesUp_uk,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goCC_diffexpGenesUp_uk<-enrichGO(gene=diffexpGenesUp_uk,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
kegg_diffexpGenesUp_uk<-enrichKEGG(diffexpGenesUp_uk_entrez,qvalueCutoff = 0.05)
reactome_diffexpGenesUp_uk<-enrichPathway(diffexpGenesUp_uk_entrez,qvalueCutoff=0.05)
enrichrWP_diffexpGenesUp_uk<-enrichr(diffexpGenesUp_uk,"WikiPathways_2019_Human")
enrichrCov2sets_diffexpGenesUp_uk<-enrichr(diffexpGenesUp_uk,"COVID-19_Related_Gene_Sets")

goBP_diffexpGenesUp_sa<-enrichGO(gene=diffexpGenesUp_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goMF_diffexpGenesUp_sa<-enrichGO(gene=diffexpGenesUp_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goCC_diffexpGenesUp_sa<-enrichGO(gene=diffexpGenesUp_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
kegg_diffexpGenesUp_sa<-enrichKEGG(diffexpGenesUp_sa_entrez,qvalueCutoff = 0.05)
reactome_diffexpGenesUp_sa<-enrichPathway(diffexpGenesUp_sa_entrez,pvalueCutoff=0.05)
enrichrWP_diffexpGenesUp_sa<-enrichr(diffexpGenesUp_sa,"WikiPathways_2019_Human")
enrichrCov2sets_diffexpGenesUp_sa<-enrichr(diffexpGenesUp_sa,"COVID-19_Related_Gene_Sets")


pdf("/plots/compiled_figures/enrichmentGO_DEGs_Upregulated_by_Strain.pdf", width = 16, height = 10)

jointdotplot(goBP_diffexpGenesUp_wu,goBP_diffexpGenesUp_uk,goBP_diffexpGenesUp_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "GO:Biological process \n Up-regulated DEGs")
jointdotplot(goMF_diffexpGenesUp_wu,goMF_diffexpGenesUp_uk,goMF_diffexpGenesUp_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "GO:Molecular function \n Up-regulated DEGs")
jointdotplot(goCC_diffexpGenesUp_wu,goCC_diffexpGenesUp_uk,goCC_diffexpGenesUp_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "GO:Cellular component \n Up-regulated DEGs")
jointdotplot(kegg_diffexpGenesUp_wu,kegg_diffexpGenesUp_uk,kegg_diffexpGenesUp_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "Kegg \n Up-regulated DEGs")

plotEnrich(enrichrCov2sets_diffexpGenesUp_wu$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))  + labs(title="Enrichr: COVID-19 Related Gene Sets \ Up-regulated genes \n B.1")
plotEnrich(enrichrCov2sets_diffexpGenesUp_uk$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))  + labs(title="Enrichr: COVID-19 Related Gene Sets \ Up-regulated genes \n B.1.1.7")
plotEnrich(enrichrCov2sets_diffexpGenesUp_sa$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))  + labs(title="Enrichr: COVID-19 Related Gene Sets \ Up-regulated genes \n B.1.351")

plotEnrich(enrichrWP_diffexpGenesUp_wu$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))   + labs(title="Enrichr: WikiPathways \ Up-regulated genes \n B.1")
plotEnrich(enrichrWP_diffexpGenesUp_uk$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))   + labs(title="Enrichr: WikiPathways \ Up-regulated genes \n B.1.1.7")
plotEnrich(enrichrWP_diffexpGenesUp_sa$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))   + labs(title="Enrichr: WikiPathways \ Up-regulated genes \n B.1.351")

dev.off()



# Supplementary Fig1C -----------------------------------------------------
pdf("/plots/compiled_figures/compilation/suppFig1C_enrichmentGO_DEGs_Up-and-Down_regulated_by_Strain.pdf", width = 16, height = 8)
plot_grid(
jointdotplot(goBP_diffexpGenesUp_wu,goBP_diffexpGenesUp_uk,goBP_diffexpGenesUp_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "GO:Biological process \n Up-regulated DEGs") + scale_y_discrete(position = "left"),
jointdotplot(goBP_diffexpGenesDown_wu,goBP_diffexpGenesDown_uk,goBP_diffexpGenesDown_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "GO:Biological process \n Down-regulated DEGs"), nrow = 1, ncol = 2
)
dev.off()

goTables_Up<-list("BP_Up_B.1"=goBP_diffexpGenesUp_wu@result %>% filter(pvalue < 0.05) ,"BP_Up_B.1.1.7"=goBP_diffexpGenesUp_uk@result %>% filter(pvalue < 0.05) ,"BP_Up_B.1.351"=goBP_diffexpGenesUp_sa@result %>% filter(pvalue < 0.05) ,"MF_Up_B.1"=goMF_diffexpGenesUp_wu@result %>% filter(pvalue < 0.05) ,"MF_Up_B.1.1.7"=goMF_diffexpGenesUp_uk@result %>% filter(pvalue < 0.05) ,"MF_Up_B.1.351"=goMF_diffexpGenesUp_sa@result %>% filter(pvalue < 0.05) ,"CC_Up_B.1"=goCC_diffexpGenesUp_wu@result %>% filter(pvalue < 0.05) ,"CC_Up_B.1.351"=goCC_diffexpGenesUp_uk@result %>% filter(pvalue < 0.05) ,"CC_Up_B.1.351"=goCC_diffexpGenesUp_sa@result %>% filter(pvalue < 0.05) ,"KEGG_Up_B.1"=kegg_diffexpGenesUp_wu@result %>% filter(pvalue < 0.05)  ,"KEGG_Up_B.1.1.7"=kegg_diffexpGenesUp_uk@result %>% filter(pvalue < 0.05) ,"KEGG_Up_B.1.351"=kegg_diffexpGenesUp_sa@result %>% filter(pvalue < 0.05) ,"EnrichrCov2Sets_B.1"=enrichrCov2sets_diffexpGenesUp_wu$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"EnrichrCov2Sets_B.1.1.7"=enrichrCov2sets_diffexpGenesUp_uk$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"EnrichrCov2Sets_B.1.351"=enrichrCov2sets_diffexpGenesUp_sa$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"KEGG_Up_B.1"=enrichrWP_diffexpGenesUp_wu$WikiPathways_2019_Human %>% filter(P.value < 0.05) ,"KEGG_Up_B.1.1.7"=enrichrWP_diffexpGenesUp_wu$WikiPathways_2019_Human %>% filter(P.value < 0.05) ,"KEGG_Up_B.1.351"=enrichrWP_diffexpGenesUp_wu$WikiPathways_2019_Human %>% filter(P.value < 0.05))


openxlsx::write.xlsx(list("Up_B.1"=diffexpGenesUp_wu,"Up_B.1.1.7"=diffexpGenesUp_uk,"Up_B.1.351"=diffexpGenesUp_sa),"/plots/compiled_figures/tables/DEGs_Up-regulated_by_Strain_geneName_list.xlsx")
openxlsx::write.xlsx(goTables_Up,"/plots/compiled_figures/tables/enrcihmenGO_DEGs_Up-regulated_by_Strain.xlsx")

#Down regulated genes:
goBP_diffexpGenesDown_wu<-enrichGO(gene=diffexpGenesDown_wu,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goMF_diffexpGenesDown_wu<-enrichGO(gene=diffexpGenesDown_wu,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goCC_diffexpGenesDown_wu<-enrichGO(gene=diffexpGenesDown_wu,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
kegg_diffexpGenesDown_wu<-enrichKEGG(diffexpGenesDown_wu_entrez,pvalueCutoff = 0.05)
reactome_diffexpGenesDown_wu<-enrichPathway(diffexpGenesDown_wu_entrez,pvalueCutoff=0.05)
enrichrWP_diffexpGenesDown_wu<-enrichr(diffexpGenesDown_wu,"WikiPathways_2019_Human")
enrichrCov2sets_diffexpGenesDown_wu<-enrichr(diffexpGenesDown_wu,"COVID-19_Related_Gene_Sets")

goBP_diffexpGenesDown_uk<-enrichGO(gene=diffexpGenesDown_uk,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goMF_diffexpGenesDown_uk<-enrichGO(gene=diffexpGenesDown_uk,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goCC_diffexpGenesDown_uk<-enrichGO(gene=diffexpGenesDown_uk,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
kegg_diffexpGenesDown_uk<-enrichKEGG(diffexpGenesDown_uk_entrez,pvalueCutoff = 0.05)
reactome_diffexpGenesDown_uk<-enrichPathway(diffexpGenesDown_uk_entrez,pvalueCutoff=0.05)
enrichrWP_diffexpGenesDown_uk<-enrichr(diffexpGenesDown_uk,"WikiPathways_2019_Human")
enrichrCov2sets_diffexpGenesDown_uk<-enrichr(diffexpGenesDown_uk,"COVID-19_Related_Gene_Sets")

goBP_diffexpGenesDown_sa<-enrichGO(gene=diffexpGenesDown_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goMF_diffexpGenesDown_sa<-enrichGO(gene=diffexpGenesDown_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goCC_diffexpGenesDown_sa<-enrichGO(gene=diffexpGenesDown_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
kegg_diffexpGenesDown_sa<-enrichKEGG(diffexpGenesDown_sa_entrez,pvalueCutoff = 0.05)
reactome_diffexpGenesDown_sa<-enrichPathway(diffexpGenesDown_sa_entrez,pvalueCutoff=0.05)
enrichrWP_diffexpGenesDown_sa<-enrichr(diffexpGenesDown_sa,"WikiPathways_2019_Human")
enrichrCov2sets_diffexpGenesDown_sa<-enrichr(diffexpGenesDown_sa,"COVID-19_Related_Gene_Sets")


pdf("/plots/compiled_figures/enrichmentGO_DEGs_Downregulated_by_Strain.pdf", width = 16, height = 10)
jointdotplot(goBP_diffexpGenesDown_wu,goBP_diffexpGenesDown_uk,goBP_diffexpGenesDown_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "GO:Biological process \n Down-regulated DEGs")
jointdotplot(goMF_diffexpGenesDown_wu,goMF_diffexpGenesDown_uk,goMF_diffexpGenesDown_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "GO:Molecular function \n Down-regulated DEGs")
jointdotplot(goCC_diffexpGenesDown_wu,goCC_diffexpGenesDown_uk,goCC_diffexpGenesDown_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "GO:Cellular component \n Down-regulated DEGs")
jointdotplot(kegg_diffexpGenesDown_wu,kegg_diffexpGenesDown_uk,kegg_diffexpGenesDown_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "Kegg \n Down-regulated DEGs")

plotEnrich(enrichrCov2sets_diffexpGenesDown_wu$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))  + labs(title="Enrichr: COVID-19 Related Gene Sets \ Down-regulated genes (B.1)")
plotEnrich(enrichrCov2sets_diffexpGenesDown_uk$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))  + labs(title="Enrichr: COVID-19 Related Gene Sets \ Down-regulated genes (B.1.1.7)")
plotEnrich(enrichrCov2sets_diffexpGenesDown_sa$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))  + labs(title="Enrichr: COVID-19 Related Gene Sets \ Down-regulated genes (B.1.351)")

plotEnrich(enrichrWP_diffexpGenesDown_wu$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))   + labs(title="Enrichr: WikiPathways \ Down-regulated genes (B.1)")
plotEnrich(enrichrWP_diffexpGenesDown_uk$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))   + labs(title="Enrichr: WikiPathways \ Down-regulated genes (B.1.1.7)")
plotEnrich(enrichrWP_diffexpGenesDown_sa$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))   + labs(title="Enrichr: WikiPathways \ Down-regulated genes (B.1.351)")

dev.off()


pdf("/plots/compiled_figures/enrichment_EnrichrCOVIDsets_UpandDown_regulated_by_Strain.pdf", width = 16, height = 10)

jointdotplotEnrichr(enrichrCov2sets_diffexpGenesDown_wu$`COVID-19_Related_Gene_Sets`,
                    enrichrCov2sets_diffexpGenesDown_uk$`COVID-19_Related_Gene_Sets`,
                    enrichrCov2sets_diffexpGenesDown_sa$`COVID-19_Related_Gene_Sets`,
                    names = c("B.1","B.1.1.7","B.1.351"),title="Enrichr: COVID-19 Down-regulated genes")

jointdotplotEnrichr(enrichrCov2sets_diffexpGenesUp_wu$`COVID-19_Related_Gene_Sets`,
                    enrichrCov2sets_diffexpGenesUp_uk$`COVID-19_Related_Gene_Sets`,
                    enrichrCov2sets_diffexpGenesUp_sa$`COVID-19_Related_Gene_Sets`,
                    names = c("B.1","B.1.1.7","B.1.351"),title="Enrichr: COVID-19 Up-regulated genes")

dev.off()

goTables_Down<-list("BP_Down_B.1"=goBP_diffexpGenesDown_wu@result %>% filter(pvalue < 0.05),"BP_Down_B.1.1.7"=goBP_diffexpGenesDown_uk@result %>% filter(pvalue < 0.05),"BP_Down_B.1.351"=goBP_diffexpGenesDown_sa@result %>% filter(pvalue < 0.05),"MF_Down_B.1"=goMF_diffexpGenesDown_wu@result %>% filter(pvalue < 0.05),"MF_Down_B.1.1.7"=goMF_diffexpGenesDown_uk@result %>% filter(pvalue < 0.05),"MF_Down_B.1.351"=goMF_diffexpGenesDown_sa@result %>% filter(pvalue < 0.05),"CC_Down_B.1"=goCC_diffexpGenesDown_wu@result %>% filter(pvalue < 0.05),"CC_Down_B.1.1.7"=goCC_diffexpGenesDown_uk@result %>% filter(pvalue < 0.05),"CC_Down_B.1.351"=goCC_diffexpGenesDown_sa@result %>% filter(pvalue < 0.05),"KEGG_Down_B.1"=kegg_diffexpGenesDown_wu@result %>% filter(pvalue < 0.05),"KEGG_Down_B.1.1.7"=kegg_diffexpGenesDown_uk@result %>% filter(pvalue < 0.05),"KEGG_Down_B.1.351"=kegg_diffexpGenesDown_sa@result %>% filter(pvalue < 0.05),"EnrichrCov2Sets_B.1"=enrichrCov2sets_diffexpGenesDown_wu$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"EnrichrCov2Sets_B.1.1.7"=enrichrCov2sets_diffexpGenesDown_uk$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"EnrichrCov2Sets_B.1.351"=enrichrCov2sets_diffexpGenesDown_sa$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"WikiPath_Down_B.1"=enrichrWP_diffexpGenesDown_wu$WikiPathways_2019_Humangenes,"WikiPath_Down_B.1.1.7"=enrichrWP_diffexpGenesDown_wu$WikiPathways_2019_Human %>% filter(P.value < 0.05),"WikiPath_Down_B.1.351"=enrichrWP_diffexpGenesDown_wu$WikiPathways_2019_Human %>% filter(P.value < 0.05))

openxlsx::write.xlsx(list("Down_B.1"=diffexpGenesDown_wu,"Down_B.1.1.7"=diffexpGenesDown_uk,"Down_B.1.351"=diffexpGenesDown_sa),"/plots/compiled_figures/tables/Down-regulatedDEGs_geneName_list.xlsx")
openxlsx::write.xlsx(goTables_Down,"/plots/compiled_figures/tables/EnrcihmenAnalysis_Down-regulatedDEGs_by_Strain.xlsx")


exclusive_wu<-setdiff(diffexpGenes_wu,diffexpGenes_uk)

#reactome_commonDiffExpGenes %>% slot("result") %>% as.data.frame() %>% arrange(pvalue,Count) %>% head(30)  %>% ggplot(aes(Count,reorder(Description,Count),size=GeneRatio,color=pvalue)) + geom_point()

# mydotplot()<-function(df,title){
#   df %>% arrange(pvalue,Count) %>%  
#   %>% ggplot(aes(x=Count,y=Description))
# }

names(diffexpGenes_wu)<-diffexpGenes_wu
names(diffexpGenes_uk)<-diffexpGenes_uk
names(diffexpGenes_sa)<-diffexpGenes_sa

ggvenn(list("B.1"=diffexpGenes_wu,"B.1.1.7"=diffexpGenes_uk,"B.1.351"=diffexpGenes_sa),columns = c("B.1","B.1.1.7","B.1.351"),stroke_color = "white",fill_color = c("#cb181d","#31a354","#e6550d")) + labs(title="Differentially expressed genes after infection") -> p1
p1
ggsave("/plots/compiled_figures/SuppFig1_vennDiagram_allDEGs_after_infection.pdf", width = 8, height = 8,dpi=300,device = "pdf", plot = p1)


names(diffexpGenesUp_wu)<-diffexpGenesUp_wu
names(diffexpGenesUp_uk)<-diffexpGenesUp_uk
names(diffexpGenesUp_sa)<-diffexpGenesUp_sa

ggvenn(list("B.1"=diffexpGenesUp_wu,"B.1.1.7"=diffexpGenesUp_uk,"B.1.351"=diffexpGenesUp_sa),columns = c("B.1","B.1.1.7","B.1.351"),stroke_color = "white",fill_color = c("#cb181d","#31a354","#e6550d")) + labs(title="Upregulated genes") -> p2
p2
ggsave("/plots/compiled_figures/vennDiagram_Up-regulatedDEGs.pdf", plot = p2 ,width = 8, height = 8,dpi=300,device = "pdf")

names(diffexpGenesDown_wu)<-diffexpGenesDown_wu
names(diffexpGenesDown_uk)<-diffexpGenesDown_uk
names(diffexpGenesDown_sa)<-diffexpGenesDown_sa

ggvenn(list("B.1"=diffexpGenesDown_wu,"B.1.1.7"=diffexpGenesDown_uk,"B.1.351"=diffexpGenesDown_sa),columns = c("B.1","B.1.1.7","B.1.351"),stroke_color = "white",fill_color = c("#cb181d","#31a354","#e6550d")) + labs(title="Downregulated genes") -> p3
p3
ggsave("/plots/compiled_figures/vennDiagram_Down-regulatedDEGs.pdf", plot = p3 ,width = 8, height = 8,dpi=300,device = "pdf")


ggsave("/plots/compiled_figures/Fig1_vennDiagram_UpandDown-regulated_DEGs.pdf", plot = plot_grid(p2,p3, labels = c("A","B"),  align = "h", vjust =7.5, label_size = 14)  ,width = 8.5, height = 8.5,dpi=300,device = "pdf")


res_all_lfcShrink %>% mutate(deg_category=if_else(log2FoldChange >1,"Up",if_else(log2FoldChange < -1,"Down",if_else(log2FoldChange<1 & log2FoldChange >-1, "NoDEG","NA")))) %>% filter(padj<0.01 &! is.na(padj) &!(gene %in% cov2genes)) %>% group_by(contrast) %>% dplyr::select(gene,deg_category) %>% distinct() %>% dplyr::count(deg_category) %>% separate(contrast,into =c("strain","non_infected"),sep="_vs_") %>% openxlsx::write.xlsx("/plots/tables/numberOf_diffexpGenes.xlsx")
#res_all_lfcShrink %>% mutate(deg_category=if_else(log2FoldChange >2,"Up",if_else(log2FoldChange < -2,"Down",if_else(log2FoldChange<2 & log2FoldChange >-2, "NoDEG","NA")))) %>% filter(padj<0.01 &! is.na(padj)) %>% group_by(contrast) %>% dplyr::select(gene,deg_category) %>% distinct() %>% dplyr::count(deg_category) %>% separate(contrast,into =c("strain","non_infected"),sep="_vs_") %>% ggplot(.,aes(strain,n,fill=deg_category)) + geom_bar(stat = "identity",position = position_stack(),width = 0.5) + scale_fill_viridis_d() + theme(panel.background = element_blank(),panel.grid = element_blank()) + theme_light()

count_diffexpGenes<-c("B.1"=length(unique(diffexpGenes_wu)),"B.1.1.7"=length(unique(diffexpGenes_uk)),"B.1.351"=length(unique(diffexpGenes_sa)))


openxlsx::write.xlsx(res_all_lfcShrink %>% filter(gene %in% common_diffexp_genes) %>% dplyr::select(gene) %>% distinct(), "/plots/tables/Table_differentialExpression_CommonDEGs_AllConditions.xls")

#normalized_counts <- counts(dds, normalized=TRUE)

commonDEGs_wu_uk <- setdiff(intersect(diffexpGenes_wu,diffexpGenes_uk),common_diffexp_genes)
commonDEGs_wu_sa <- setdiff(intersect(diffexpGenes_wu,diffexpGenes_sa),common_diffexp_genes)
commonDEGs_uk_sa <- setdiff(intersect(diffexpGenes_uk,diffexpGenes_sa),common_diffexp_genes)
names(commonDEGs_wu_uk)<-commonDEGs_wu_uk
names(commonDEGs_wu_sa)<-commonDEGs_wu_sa
names(commonDEGs_uk_sa)<-commonDEGs_uk_sa

commonDEGs_wu_uk_entrez<-mapIds(org.Hs.eg.db,keys = commonDEGs_wu_uk, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
commonDEGs_wu_sa_entrez<-mapIds(org.Hs.eg.db,keys = commonDEGs_wu_sa, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
commonDEGs_uk_sa_entrez<-mapIds(org.Hs.eg.db,keys = commonDEGs_uk_sa, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
names(commonDEGs_wu_uk_entrez)<-commonDEGs_wu_uk_entrez
names(commonDEGs_wu_sa_entrez)<-commonDEGs_wu_sa_entrez
names(commonDEGs_uk_sa_entrez)<-commonDEGs_uk_sa_entrez



goBP_commonDEGs_wu_uk<-enrichGO(gene=commonDEGs_wu_uk,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goMF_commonDEGs_wu_uk<-enrichGO(gene=commonDEGs_wu_uk,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goCC_commonDEGs_wu_uk<-enrichGO(gene=commonDEGs_wu_uk,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
kegg_commonDEGs_wu_uk<-enrichKEGG(commonDEGs_wu_uk_entrez,pvalueCutoff = 0.05)
reactome_commonDEGs_wu_uk<-enrichPathway(commonDEGs_wu_uk_entrez,pvalueCutoff=0.05)
enrichrWP_commonDEGs_wu_uk<-enrichr(commonDEGs_wu_uk,"WikiPathways_2019_Human")
enrichrCov2sets_commonDEGs_wu_uk<-enrichr(commonDEGs_wu_uk,"COVID-19_Related_Gene_Sets")

goBP_commonDEGs_wu_sa<-enrichGO(gene=commonDEGs_wu_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goMF_commonDEGs_wu_sa<-enrichGO(gene=commonDEGs_wu_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goCC_commonDEGs_wu_sa<-enrichGO(gene=commonDEGs_wu_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
kegg_commonDEGs_wu_sa<-enrichKEGG(commonDEGs_wu_sa_entrez,pvalueCutoff = 0.05)
reactome_commonDEGs_wu_sa<-enrichPathway(commonDEGs_wu_sa_entrez,pvalueCutoff=0.05)
enrichrWP_commonDEGs_wu_sa<-enrichr(commonDEGs_wu_sa,"WikiPathways_2019_Human")
enrichrCov2sets_commonDEGs_wu_sa<-enrichr(commonDEGs_wu_sa,"COVID-19_Related_Gene_Sets")

goBP_commonDEGs_uk_sa<-enrichGO(gene=commonDEGs_uk_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goMF_commonDEGs_uk_sa<-enrichGO(gene=commonDEGs_uk_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goCC_commonDEGs_uk_sa<-enrichGO(gene=commonDEGs_uk_sa,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
kegg_commonDEGs_uk_sa<-enrichKEGG(commonDEGs_uk_sa_entrez,pvalueCutoff = 0.05)
reactome_commonDEGs_uk_sa<-enrichPathway(commonDEGs_uk_sa_entrez,pvalueCutoff=0.05)
enrichrWP_commonDEGs_uk_sa<-enrichr(commonDEGs_uk_sa,"WikiPathways_2019_Human")
enrichrCov2sets_commonDEGs_uk_sa<-enrichr(commonDEGs_uk_sa,"COVID-19_Related_Gene_Sets")


pdf("/plots/compiled_figures/enrichmentGO_commonDEGs_Paired_Strain_Comparisons.pdf", width = 20, height = 10)

jointdotplot(goBP_commonDEGs_wu_uk,goBP_commonDEGs_wu_sa,goBP_commonDEGs_uk_sa,names = c("B.1\n&\nB.1.1.7","B.1 & B.1.351","B.1.1.7 & B.1.351"),title = "GO:Biological Process \n Common DEGs Paired Strain Comparisons")
jointdotplot(goMF_commonDEGs_wu_uk,goMF_commonDEGs_wu_sa,goMF_commonDEGs_uk_sa,names = c("B.1 & B.1.1.7","B.1 & B.1.351","B.1.1.7 & B.1.351"),title = "GO:Molecular Function \n Common DEGs Paired Strain Comparisons")
jointdotplot(goCC_commonDEGs_wu_uk,goCC_commonDEGs_wu_sa,goCC_commonDEGs_uk_sa,names = c("B.1 & B.1.1.7","B.1 & B.1.351","B.1.1.7 & B.1.351"),title = "GO:Cellular Component \n Common DEGs Paired Strain Comparisons")
jointdotplot(kegg_commonDEGs_wu_uk,kegg_commonDEGs_wu_sa,kegg_commonDEGs_uk_sa,names = c("B.1 & B.1.1.7","B.1 & B.1.351","B.1.1.7 & B.1.351"),title = "KEGG \n Common DEGs Paired Strain Comparisons")

plotEnrich(enrichrCov2sets_commonDEGs_wu_uk$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))  + labs(title="Enrichr: COVID-19 Related Gene Sets \ Common DEGs B.1 & B.1.1.7)")
plotEnrich(enrichrCov2sets_commonDEGs_wu_sa$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))  + labs(title="Enrichr: COVID-19 Related Gene Sets \ Common DEGs (B.1 & B.1.351)")
plotEnrich(enrichrCov2sets_commonDEGs_uk_sa$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))  + labs(title="Enrichr: COVID-19 Related Gene Sets \ Common DEGs (B.1.351 & B.1.1.7)")

plotEnrich(enrichrWP_commonDEGs_wu_uk$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))  + labs(title="Enrichr: WikiPathways \ Common DEGs (B.1 & B.1.1.7)")
plotEnrich(enrichrWP_commonDEGs_wu_sa$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))  + labs(title="Enrichr: WikiPathways \ Common DEGs (B.1 & B.1.351)")
plotEnrich(enrichrWP_commonDEGs_uk_sa$WikiPathways_2019_Human %>% filter(P.value<0.05),numChar = 100) + theme(aspect.ratio = 2.5,axis.text.y = element_text(size=7))  + labs(title="Enrichr: WikiPathways \ Common DEGs (B.1.351 & B.1.1.7)")

dev.off()

goTables_commonDEgs_pairedStrains<-list("BP_commonDEGs_WU_UK"=goBP_commonDEGs_wu_uk@result %>% filter(pvalue < 0.05),"BP_commonDEGs_SA_UK"=goBP_commonDEGs_uk_sa@result %>% filter(pvalue < 0.05),"BP_commonDEGs_WU_SA"=goBP_commonDEGs_wu_sa@result %>% filter(pvalue < 0.05),"MF_commonDEGs_WU_UK"=goMF_commonDEGs_wu_uk@result %>% filter(pvalue < 0.05),"MF_commonDEGs_SA_UK"=goMF_commonDEGs_uk_sa@result %>% filter(pvalue < 0.05),"MF_commonDEGs_WU_SA"=goMF_commonDEGs_wu_sa@result %>% filter(pvalue < 0.05),"CC_commonDEGs_WU_UK"=goCC_commonDEGs_wu_uk@result %>% filter(pvalue < 0.05),"CC_commonDEGs_SA_UK"=goCC_commonDEGs_uk_sa@result %>% filter(pvalue < 0.05),"CC_commonDEGs_WU_SA"=goCC_commonDEGs_wu_sa@result %>% filter(pvalue < 0.05),"KEGG_commonDEGs_WU_UK"=kegg_commonDEGs_wu_uk@result %>% filter(pvalue < 0.05),"KEGG_commonDEGs_SA_UK"=kegg_commonDEGs_uk_sa@result %>% filter(pvalue < 0.05),"KEGG_commonDEGs_WU_SA"=kegg_commonDEGs_wu_sa@result %>% filter(pvalue < 0.05),"EnrichrCov2Sets_WU"=enrichrCov2sets_commonDEGs_wu_uk$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"EnrichrCov2Sets_UK"=enrichrCov2sets_commonDEGs_uk_sa$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"EnrichrCov2Sets_SA"=enrichrCov2sets_commonDEGs_wu_sa$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),"WikiPath_commonDEGs_WU_UK"=enrichrWP_commonDEGs_wu_uk$WikiPathways_2019_Humangenes,"WikiPath_commonDEGs_SA_UK"=enrichrWP_commonDEGs_wu_uk$WikiPathways_2019_Human %>% filter(P.value < 0.05),"WikiPath_commonDEGs_WU_SA"=enrichrWP_commonDEGs_wu_uk$WikiPathways_2019_Human %>% filter(P.value < 0.05))

openxlsx::write.xlsx(list("Down_WU"=commonDEGs_wu_uk,"Down_UK"=commonDEGs_uk_sa,"Down_SA"=commonDEGs_wu_sa),"/plots/tables/commonDEGs_Paired_Strain_Comparisons_list.xlsx")
openxlsx::write.xlsx(goTables_commonDEgs_pairedStrains,"/plots/tables/EnrcihmenAnalysis_commonDEGs_Paired_Strain_Comparisons.xlsx")



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


#library(ggtext)

jointdotplot<-function(df1,df2,df3,names,title){

  #file<-tolower(gsub(" ","",file))
  
  df1@result$GeneRatio<-sapply(df1@result$GeneRatio,function(x){eval(parse(text=x))})
  df2@result$GeneRatio<-sapply(df2@result$GeneRatio,function(x){eval(parse(text=x))})
  df3@result$GeneRatio<-sapply(df3@result$GeneRatio,function(x){eval(parse(text=x))})
  
  df1_filt<-df1@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[1], compose_score=-log10(pvalue)*GeneRatio) %>%  arrange(desc(compose_score)) %>% head(30)
  df2_filt<-df2@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[2],compose_score=-log10(pvalue)*GeneRatio) %>%  arrange(desc(compose_score)) %>% head(30)
  df3_filt<-df3@result %>% filter(pvalue < 0.05 ) %>% mutate(condition=names[3], compose_score=-log10(pvalue)*GeneRatio) %>%  arrange(desc(compose_score))%>% head(30)
  
  topTerms<-unique(union(union(df1_filt$Description,df2_filt$Description),df3_filt$Description))
  commonTerms<-unique(intersect(intersect(df1_filt$Description,df2_filt$Description),df3_filt$Description))
  diffTerms<-setdiff(topTerms,commonTerms)
  
  #df1_filt<-df1@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[1])
  #df2_filt<-df2@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[2])
  #df3_filt<-df3@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[3])
  
  df<-rbind(df1_filt,df2_filt,df3_filt)
  #df$condition <- factor(df$condition,levels = names)
  df_filt<- df %>% filter(Description %in% topTerms) %>% arrange(desc(-log10(pvalue)),desc(GeneRatio))
  df_wide <-df_filt %>% pivot_wider(names_from = "condition",values_from = colnames(df_filt)[-c(1,2,10)])
  
  #diffTerms <- df_wide %>% mutate(is_diff=if_else((is.na(GeneRatio_WU) | is.na(GeneRatio_UK) |is.na(GeneRatio_SA)),TRUE,FALSE)) %>% filter(is_diff==TRUE) %>% dplyr::select(Description) %>% distinct() %>% unlist()
  # df_filt<- left_join(df_filt,df_filt %>% count(Description),by="Description") %>%  
  #   mutate(labelcolor=if_else(n==3,glue("<span style='color:#000000'>{Description}</span>"),
  #                             glue("<span style='color:#d7301f'>{Description}</span>"),
  #                             glue("<span style='color:#000000;'>{Description}</span>")))
  df_filt<- df_filt %>% arrange(dplyr::desc(compose_score))
  #df_filt$Description
  #df_filt$Description<-factor(df_filt$Description,levels = unique(df_filt$Description))
  #order<-df_filt$Description
  
  #message("Saving top terms to file")
  #openxlsx::write.xlsx(x = df_wide,file = paste0(outdir,"/",file))
  #aes(fct_relevel(condition,levels=c("B.1","B.1.1.7","B.1.351")),reorder(Description,-log10(pvalue)*GeneRatio)
  
   df_filt %>% ggplot(.,aes(fct_relevel(condition,levels=c("B.1","B.1.1.7","B.1.351")),reorder(Description,compose_score), size=GeneRatio,color=-log10(pvalue))) + 
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


jointdotplot2<-function(df1,df2,df3,names,title){
  
  #file<-tolower(gsub(" ","",file))
  
  df1@result$GeneRatio<-sapply(df1@result$GeneRatio,function(x){eval(parse(text=x))})
  df2@result$GeneRatio<-sapply(df2@result$GeneRatio,function(x){eval(parse(text=x))})
  df3@result$GeneRatio<-sapply(df3@result$GeneRatio,function(x){eval(parse(text=x))})
  
  df1_filt<-df1@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[1], compose_score=-log10(pvalue)*GeneRatio) %>%  arrange(desc(compose_score)) %>% head(30)
  df2_filt<-df2@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[2],compose_score=-log10(pvalue)*GeneRatio) %>%  arrange(desc(compose_score)) %>% head(30)
  df3_filt<-df3@result %>% filter(pvalue < 0.05 ) %>% mutate(condition=names[3], compose_score=-log10(pvalue)*GeneRatio) %>%  arrange(desc(compose_score))%>% head(30)
  
  topTerms<-unique(union(union(df1_filt$Description,df2_filt$Description),df3_filt$Description))
  commonTerms<-unique(intersect(intersect(df1_filt$Description,df2_filt$Description),df3_filt$Description))
  diffTerms<-setdiff(topTerms,commonTerms)
  
  #df1_filt<-df1@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[1])
  #df2_filt<-df2@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[2])
  #df3_filt<-df3@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[3])
  
  df<-rbind(df1_filt,df2_filt,df3_filt)
  #df$condition <- factor(df$condition,levels = names)
  df_filt<- df %>% filter(Description %in% topTerms) %>% arrange(desc(-log10(pvalue)),desc(GeneRatio))
  df_wide <-df_filt %>% pivot_wider(names_from = "condition",values_from = colnames(df_filt)[-c(1,2,10)])
  
  #diffTerms <- df_wide %>% mutate(is_diff=if_else((is.na(GeneRatio_WU) | is.na(GeneRatio_UK) |is.na(GeneRatio_SA)),TRUE,FALSE)) %>% filter(is_diff==TRUE) %>% dplyr::select(Description) %>% distinct() %>% unlist()
  # df_filt<- left_join(df_filt,df_filt %>% count(Description),by="Description") %>%  
  #   mutate(labelcolor=if_else(n==3,glue("<span style='color:#000000'>{Description}</span>"),
  #                             glue("<span style='color:#d7301f'>{Description}</span>"),
  #                             glue("<span style='color:#000000;'>{Description}</span>")))
  df_filt<- df_filt %>% arrange(dplyr::desc(compose_score))
  #df_filt$Description
  #df_filt$Description<-factor(df_filt$Description,levels = unique(df_filt$Description))
  #order<-df_filt$Description
  
  #message("Saving top terms to file")
  #openxlsx::write.xlsx(x = df_wide,file = paste0(outdir,"/",file))
  #aes(fct_relevel(condition,levels=c("B.1","B.1.1.7","B.1.351")),reorder(Description,-log10(pvalue)*GeneRatio)
  
  df_filt %>% ggplot(.,aes(fct_relevel(condition,levels=c("B.1","B.1.1.7","B.1.351")),reorder(Description,compose_score), size=GeneRatio,color=-log10(pvalue))) + 
    geom_point(alpha=0.8) +
    scale_y_discrete(position = "right") +
    theme_light() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_markdown(angle=0,hjust = 0.5,size=8),
          axis.text.y = element_markdown(size=6),aspect.ratio = 5) +
    scale_color_gradient(high="salmon",low = "blue") +
    labs(title=title) + xlab("") + ylab("") +
    coord_flip()->p
  return(p)
}

library(ggtext)
#Heatmap
library(pheatmap)
res_diffexp_df<- res_all_lfcShrink %>% filter(abs(log2FoldChange)>1 & padj < 0.05 & !is.na(padj))
res_diffexp_df <- res_diffexp_df %>% dplyr::select(symbol,log2FoldChange,contrast) %>%
  pivot_wider(names_from = contrast,values_from = log2FoldChange) %>%  
  mutate(wu_vs_vero=replace_na(wu_vs_vero,replace = 0),uk_vs_vero=replace_na(uk_vs_vero,replace = 0),sa_vs_vero=replace_na(sa_vs_vero,replace = 0))

res_diffexp.m<-as.matrix(res_diffexp_df[,-1])
rownames(res_diffexp.m)<-res_diffexp_df[,1]

pheatmap(res_diffexp.m,cluster_rows = T,clus)



# Common GO terms BP for ClueGo analysis ----------------------------------
common_goBP<-inner_join(go_bp_sa@result %>% filter(pvalue<0.05),
                        inner_join(go_bp_wu@result %>% filter(pvalue<0.05),go_bp_uk@result %>% filter(pvalue<0.05),by="ID"),by="ID") %>%
  dplyr::select(ID,pvalue,Description)
union_goBP<-rbind(go_bp_wu@result %>% filter(pvalue<0.05) %>% dplyr::select(ID,pvalue,Description, geneID),
                  go_bp_uk@result %>% filter(pvalue<0.05) %>% dplyr::select(ID,pvalue,Description,geneID),
                  go_bp_sa@result %>% filter(pvalue<0.05) %>% dplyr::select(ID,pvalue,Description, geneID)) %>% arrange(pvalue)
  

union_goBP %>%
  filter(grepl("METTL3|METTL14|KIAA1429|CBLL1|WTAP|RBM15|RBM15B|ZC3H13|SPEN||YTHDC1|YTHDC2|YTHDF1|YTHDF2|YTHDF3|HNRNPA2B1|HNRNPC|FTO|ALKBH5",geneID, perl = T)) %>% 
  filter(grepl("(viral entry into host cell | virus | viral|splic|SRP|RNA catabolic process|RNA stabilization|stress granule|translation | protein targeting | cilium | localization to membrane|translational initiation|stress response|response to stress| nuclear export|export from nucleus)",Description,perl=T,ignore.case = T)) %>%
  arrange(pvalue) %>% dplyr::select(ID,pvalue,Description) %>%
  fwrite(.,"/plots/tables/GOTerms_writers_readers_erasers_GOBP_filtered.txt", sep="\t",col.names = F)


diffExp_m6aWritersReadersErasers <- res_all_lfcShrink %>% filter(gene %in% c(writers,readers,erasers)) %>% filter(abs(log2FoldChange)>1 & padj<0.01 &!is.na(padj)) %>% dplyr::select(gene) %>% distinct() %>% unlist()

#Select enriched GO terms shared across any two strains combination:
combinations<-combn(c("B.1","B.1.1.7","B.1.351"),2)


sigTerms_wu <- gofilter(go_bp_wu,level = c(7:10))
sigTerms_uk <- gofilter(go_bp_uk,level = c(7:10))
sigTerms_sa <- gofilter(go_bp_sa,level = c(7:10))

sigTerms_wu <- go_bp_wu@result %>% filter(pvalue<0.05)
colnames(sigTerms_wu) <- paste0(colnames(sigTerms_wu),"_wu")

sigTerms_uk <- go_bp_uk@result %>% filter(pvalue<0.05)
colnames(sigTerms_uk) <- paste0(colnames(sigTerms_uk),"_uk")

sigTerms_sa <- go_bp_sa@result %>% filter(pvalue<0.05)
colnames(sigTerms_sa) <- paste0(colnames(sigTerms_sa),"_sa")

commonTerms_WU_UK <- inner_join(sigTerms_wu,sigTerms_uk, by=c("ID_wu"="ID_uk"))

commonTerms_WU_SA <- inner_join(sigTerms_wu,sigTerms_sa, by=c("ID_wu"="ID_sa"))

commonTerms_SA_UK <- inner_join(sigTerms_sa,sigTerms_uk, by=c("ID_sa"="ID_uk"))



commonTerms_WU_UK_m6Aprots <- commonTerms_WU_UK %>% 
  filter(grepl("(HNRNPC|HNRNPA2B1|SPEN|RBM15|WTAP|KIAA1429|YTHDC1|YTHDF1|FTO|RBM15B)",geneID_wu, perl=T) & grepl("(HNRNPC|HNRNPA2B1|SPEN|RBM15|WTAP|KIAA1429|YTHDC1|YTHDF1|FTO|RBM15B)",geneID_uk, perl=T))

commonTerms_WU_SA_m6Aprots <- commonTerms_WU_SA %>% 
  filter(grepl("(HNRNPC|HNRNPA2B1|SPEN|RBM15|WTAP|KIAA1429|YTHDC1|YTHDF1|FTO|RBM15B)",geneID_wu, perl=T) & grepl("(HNRNPC|HNRNPA2B1|SPEN|RBM15|WTAP|KIAA1429|YTHDC1|YTHDF1|FTO|RBM15B)",geneID_sa, perl=T))

commonTerms_SA_UK_m6Aprots <- commonTerms_SA_UK %>% 
  filter(grepl("(HNRNPC|HNRNPA2B1|SPEN|RBM15|WTAP|KIAA1429|YTHDC1|YTHDF1|FTO|RBM15B)",geneID_sa, perl=T) & grepl("(HNRNPC|HNRNPA2B1|SPEN|RBM15|WTAP|KIAA1429|YTHDC1|YTHDF1|FTO|RBM15B)",geneID_uk, perl=T))

openxlsx::write.xlsx(list("commonTerms_WU_UK"=commonTerms_WU_UK,"commonTerms_WU_SA"=commonTerms_WU_SA,"commonTerms_SA_UK"=commonTerms_SA_UK), file = "/plots/tables/commonGOTerms_biological_process_strain_pairs_All.xlsx")

openxlsx::write.xlsx(list("commonTerms_WU_UK"=commonTerms_WU_UK_m6Aprots,"commonTerms_WU_SA"=commonTerms_WU_SA_m6Aprots,"commonTerms_SA_UK"=commonTerms_SA_UK_m6Aprots), file = "/plots/tables/commonGOTerms_biological_process_strain_pairs_m6AProteinsInBothStrains.xlsx")

commonTerms_All<- inner_join(inner_join(sigTerms_wu,sigTerms_uk, by=c("ID_wu"="ID_uk")),sigTerms_sa,by=c("ID_wu"="ID_sa")) 

commonTerms_All_m6aprots <- commonTerms_All %>% 
  filter(grepl("(HNRNPC|HNRNPA2B1|SPEN|RBM15|WTAP|KIAA1429|YTHDC1|YTHDF1|FTO|RBM15B)",geneID_wu, perl=T) & grepl("(HNRNPC|HNRNPA2B1|SPEN|RBM15|WTAP|KIAA1429|YTHDC1|YTHDF1|FTO|RBM15B)",geneID_sa, perl=T))

jointdotplotHoriz<-function(df1,df2,df3,names,title){
  
  #file<-tolower(gsub(" ","",file))
  
  df1@result$GeneRatio<-sapply(df1@result$GeneRatio,function(x){eval(parse(text=x))})
  df2@result$GeneRatio<-sapply(df2@result$GeneRatio,function(x){eval(parse(text=x))})
  df3@result$GeneRatio<-sapply(df3@result$GeneRatio,function(x){eval(parse(text=x))})
  
  df1_filt<-df1@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[1], compose_score=-log10(pvalue)*GeneRatio) %>%  arrange(desc(compose_score)) %>% head(30)
  df2_filt<-df2@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[2],compose_score=-log10(pvalue)*GeneRatio) %>%  arrange(desc(compose_score)) %>% head(30)
  df3_filt<-df3@result %>% filter(pvalue < 0.05 ) %>% mutate(condition=names[3], compose_score=-log10(pvalue)*GeneRatio) %>%  arrange(desc(compose_score))%>% head(30)
  
  topTerms<-unique(union(union(df1_filt$Description,df2_filt$Description),df3_filt$Description))
  commonTerms<-unique(intersect(intersect(df1_filt$Description,df2_filt$Description),df3_filt$Description))
  diffTerms<-setdiff(topTerms,commonTerms)
  
  #df1_filt<-df1@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[1])
  #df2_filt<-df2@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[2])
  #df3_filt<-df3@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[3])
  
  df<-rbind(df1_filt,df2_filt,df3_filt)
  #df$condition <- factor(df$condition,levels = names)
  df_filt<- df %>% filter(Description %in% topTerms) %>% arrange(desc(-log10(pvalue)),desc(GeneRatio))
  df_wide <-df_filt %>% pivot_wider(names_from = "condition",values_from = colnames(df_filt)[-c(1,2,10)])
  
  #diffTerms <- df_wide %>% mutate(is_diff=if_else((is.na(GeneRatio_WU) | is.na(GeneRatio_UK) |is.na(GeneRatio_SA)),TRUE,FALSE)) %>% filter(is_diff==TRUE) %>% dplyr::select(Description) %>% distinct() %>% unlist()
  # df_filt<- left_join(df_filt,df_filt %>% count(Description),by="Description") %>%  
  #   mutate(labelcolor=if_else(n==3,glue("<span style='color:#000000'>{Description}</span>"),
  #                             glue("<span style='color:#d7301f'>{Description}</span>"),
  #                             glue("<span style='color:#000000;'>{Description}</span>")))
  df_filt<- df_filt %>% arrange(dplyr::desc(compose_score))
  #df_filt$Description
  #df_filt$Description<-factor(df_filt$Description,levels = unique(df_filt$Description))
  #order<-df_filt$Description
  
  #message("Saving top terms to file")
  #openxlsx::write.xlsx(x = df_wide,file = paste0(outdir,"/",file))
  #aes(fct_relevel(condition,levels=c("B.1","B.1.1.7","B.1.351")),reorder(Description,-log10(pvalue)*GeneRatio)
  
  df_filt %>% ggplot(.,aes(fct_relevel(condition,levels=rev(c("B.1","B.1.1.7","B.1.351"))),reorder(Description,compose_score), size=GeneRatio,color=-log10(pvalue))) + 
    geom_point(alpha=0.8) +
    scale_y_discrete(position = "left", expand = c(0.25,0.25)) +
    scale_x_discrete(position = "top") +
    theme_light() +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.text.x = element_markdown(angle=45,hjust = 1,size=8),
          axis.text.y = element_markdown(size=8),aspect.ratio = 0.1) +
    scale_color_gradient(high="salmon",low = "blue") +
    labs(title=title) + xlab("") + ylab("") +
    coord_flip() ->p
  return(p)
}

pdf("/plots/compiled_figures/compilation/Fig1C_GOBiologicalProcess_DEGsAll_horizontal_20211004.pdf", width = 12,height = 6)

jointdotplotHoriz(go_bp_wu,go_bp_uk,go_bp_sa,names = c("B.1","B.1.1.7","B.1.351"),title = "GO:Biological process \n DEGs B.1, B.1.1.7, B.1.351")

dev.off()



jointdotplotEnrichr<-function(df1,df2,df3,names,title){
  
  #file<-tolower(gsub(" ","",file))
  
  #df1@result$GeneRatio<-sapply(df1@result$GeneRatio,function(x){eval(parse(text=x))})
  #df2@result$GeneRatio<-sapply(df2@result$GeneRatio,function(x){eval(parse(text=x))})
  #df3@result$GeneRatio<-sapply(df3@result$GeneRatio,function(x){eval(parse(text=x))})
  
  df1_filt<-df1 %>% filter(P.value < 0.05 )  %>% mutate(condition=names[1], compose_score= Combined.Score) %>%  arrange(desc(compose_score)) %>% head(30)
  df2_filt<-df2 %>% filter(P.value < 0.05 )  %>% mutate(condition=names[2],compose_score= Combined.Score) %>%  arrange(desc(compose_score)) %>% head(30)
  df3_filt<-df3 %>% filter(P.value < 0.05 ) %>% mutate(condition=names[3], compose_score= Combined.Score) %>%  arrange(desc(compose_score)) %>% head(30)
  %>% head(60)
  
  topTerms<-unique(union(union(df1_filt$Term,df2_filt$Term),df3_filt$Term))
  commonTerms<-unique(intersect(intersect(df1_filt$Term,df2_filt$Term),df3_filt$Term))
  diffTerms<-setdiff(topTerms,commonTerms)
  
  #df1_filt<-df1@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[1])
  #df2_filt<-df2@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[2])
  #df3_filt<-df3@result %>% filter(pvalue < 0.05 )  %>% mutate(condition=names[3])
  
  df<-rbind(df1_filt,df2_filt,df3_filt)
  #df$condition <- factor(df$condition,levels = names)
  df_filt<- df %>% filter(Term %in% topTerms) %>% arrange(desc(-log10(P.value)),desc(Odds.Ratio))
  df_wide <-df_filt %>% pivot_wider(names_from = "condition",values_from = colnames(df_filt)[-c(1,2,10)])
  
  #diffTerms <- df_wide %>% mutate(is_diff=if_else((is.na(GeneRatio_WU) | is.na(GeneRatio_UK) |is.na(GeneRatio_SA)),TRUE,FALSE)) %>% filter(is_diff==TRUE) %>% dplyr::select(Description) %>% distinct() %>% unlist()
  # df_filt<- left_join(df_filt,df_filt %>% count(Description),by="Description") %>%  
  #   mutate(labelcolor=if_else(n==3,glue("<span style='color:#000000'>{Description}</span>"),
  #                             glue("<span style='color:#d7301f'>{Description}</span>"),
  #                             glue("<span style='color:#000000;'>{Description}</span>")))
  df_filt<- df_filt %>% arrange(dplyr::desc(compose_score))
  #df_filt$Description
  #df_filt$Description<-factor(df_filt$Description,levels = unique(df_filt$Description))
  #order<-df_filt$Description
  
  #message("Saving top terms to file")
  #openxlsx::write.xlsx(x = df_wide,file = paste0(outdir,"/",file))
  #aes(fct_relevel(condition,levels=c("B.1","B.1.1.7","B.1.351")),reorder(Description,-log10(pvalue)*GeneRatio)
  
  df_filt %>% ggplot(.,aes(fct_relevel(condition,levels=c("B.1","B.1.1.7","B.1.351")),reorder(Term,compose_score), size=Odds.Ratio,color=-log10(P.value))) + 
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

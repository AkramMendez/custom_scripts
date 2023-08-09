library(data.table)
library(tidyverse)
library(tximport)
library(rnaseqDTU)
library(GenomicFeatures)
library(DRIMSeq)
library(stageR)
library(DEXSeq)
library(cowplot)


# DEXseq data preparation -------------------------------------------------


indir<-"/dexseq_counts"
countFiles<-list.files(indir,pattern = ".txt$",full.names = T)

flattenedFile= list.files(indir, pattern = "ensembl.gff$",full.names = T)
basename(flattenedFile)
sampleTable = data.frame(
  row.names = c("SA_i1","SA_i2","UK_i1","UK_i2","Vero_i1","Vero_i2","WU_i1","WU_i2"),
  condition= c("SA","SA","UK","UK","Vero","Vero","WU","WU"),
  libType=rep("single-end",8)
)

sampleTable$condition<-factor(sampleTable$condition)
sampleTable$condition<-relevel(sampleTable$condition,"Vero")

dxd_dataset <- DEXSeqDataSetFromHTSeq(countfiles = countFiles,
                                      sampleData = sampleTable,
                                      design = ~sample + exon + condition:exon,
                                      flattenedfile = flattenedFile)


# DEXseq Normalization ----------------------------------------------------


dxd_dataset = estimateSizeFactors(dxd_dataset,locfunc=shorth)
dxd_dataset = estimateDispersions(dxd_dataset)
#par(mar=c(1,1,1,1))
plotDispEsts(dxd_dataset)



# Test for differential exon usage: ---------------------------------------


dxd_dataset = testForDEU(dxd_dataset)

dxd_dataset = estimateExonFoldChanges(dxd_dataset, fitExpToVar = "condition")


dex_results <- DEXSeqResults(dxd_dataset)
save(dex_results,file = "/dexseq_counts/dex_results_by_strain.rda")
load(file = "/dexseq_counts/dex_results_by_strain.rda")
mcols(dex_results)$description

q.cutoff<-0.1
#How many exonic regions are significant with a FDR < q.cutoff ?:
table ( dex_results$padj < q.cutoff)
#FALSE  TRUE 
#31490   954

#How many genes have diff. exon usage?:
#length(unique(names(perGeneQValue(dex_results))[perGeneQValue(dex_results)< q.cutoff]))
#740
table ( tapply( dex_results$padj < q.cutoff, dex_results$groupID, any ) )

#Controlling FDR at gene level, how many genes pass the FDR threshold?:
numOfGenes<-sum(perGeneQValue(dex_results) < q.cutoff)
numOfGenes

#Visualization
DEXSeq::plotMA(dex_results,cex=0.75) 


plotDEXSeq(dex_results,"ENSCSAG00000007240",legend = TRUE,cex.axis=1.2,cex=1.3,lwd=2, displayTranscripts=F, norCounts=FALSE, splicing=TRUE,expression = T)


pdf("ENO3_dexSEqplot.pdf")
plotDEXSeq(dex_results,"ENSCSAG00000012928",legend = TRUE,cex.axis=1.2,cex=1.3,lwd=2, displayTranscripts=TRUE, norCounts=FALSE, splicing=TRUE,expression = FALSE)
dev.off()

fwrite(dex_results %>% as.data.frame() %>% filter(padj < q.cutoff),"/tables/diffExonUsage_DEXseq_results_cov2_padj0.05.tsv", sep="\t")

dex_filt<- dex_results %>% as.data.frame() %>% filter(padj < q.cutoff)
#Convert Ensembl IDs to gene symbols:
library(biomaRt)
#httr::set_config(httr::config(ssl_verifypeer = FALSE))
mart <- useDataset("csabaeus_gene_ensembl", useMart("ensembl"))
#Fetch gene symbols from Biomart:
dex_filt_ensembl <- gsub("\\+.*","",dex_filt$groupID,perl = T)
dex_filt_gene_symbol <- getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol","entrezgene_accession","uniprot_gn_symbol","external_gene_name"),values = dex_filt_ensembl,mart = mart)
#Impute empty gene names:

dex_filt_gene_symbol<- dex_filt_gene_symbol %>% mutate(gene=if_else(hgnc_symbol=="",uniprot_gn_symbol,hgnc_symbol)) %>% mutate(gene=if_else(gene=="",ensembl_gene_id,gene,ensembl_gene_id)) 

#Join imputed gene names with DEX results data frame:
dex_filt <- left_join(dex_filt %>% mutate(groupID=gsub("\\+.*","",groupID,perl = T)),dex_filt_gene_symbol %>% dplyr::select(ensembl_gene_id,gene), by=c("groupID"="ensembl_gene_id"))
save(dex_filt, file = "/dexseq_counts/dex_filt_deuGenes_dataset_ensemblIDs_geneIDs.rda")

# Summary, number of included and included exons by condition -------------


dex_filt_includedExcluded<- dex_filt %>%
  mutate(name=paste0(dex_filt$gene,"_",dex_filt$featureID),
         WU=if_else(log2fold_WU_Vero > 0,1,0),
         UK=if_else(log2fold_UK_Vero > 0,1,0),
         SA=if_else(log2fold_SA_Vero > 0,1,0)) %>%
  dplyr::select(name,WU,UK,SA) %>% distinct()

summary_DEU_includedExcluded<- dex_filt_includedExcluded %>%
  pivot_longer(names_to = "condition", values_to = "diff_exon_usage",-name) %>%
  mutate(diff_exon_usage=if_else(diff_exon_usage==1,"included",if_else(diff_exon_usage==0,"excluded","NA")))

summary_DEU_includedExcluded %>% group_by(condition) %>% 
  dplyr::count(diff_exon_usage) %>% 
  pivot_wider(names_from = "diff_exon_usage",values_from = "n") %>%
  kbl() %>% kable_paper(full_width=F)

dex_filt_includedExcluded.m <- dex_filt_includedExcluded %>% 
  dplyr::select(WU,UK,SA) %>% as.matrix()

rownames(dex_filt_includedExcluded.m)<-dex_filt_includedExcluded$name

# Heatmap DEU exons by inclusion/exclusion status -------------------------

pdf("heatmap_includedExcluded_exons_DEU.pdf")

pheatmap(dex_filt_includedExcluded.m, fontsize_row = 1,fontsize_col = 3, cluster_rows = T, cluster_cols = F)

dev.off()


# Heatmap DEU exons by logFC values ---------------------------------------


dex_filt.m <-dex_filt %>% dplyr::select(log2fold_WU_Vero,log2fold_UK_Vero,log2fold_SA_Vero) %>% as.matrix()

rownames(dex_filt.m)<-paste0(dex_filt$gene,"_",dex_filt$featureID)

topVarGenes <- head(order(-rowVars(dex_filt.m)),30)

mat<- dex_filt.m[topVarGenes,]
mat<- mat - rowMeans(mat)
rownames(mat)<-rownames(dex_filt.m)[topVarGenes]

pdf(file="heatmap_DiffExonUsage.pdf")

pheatmap(mat)
#dex_filt$gene_symbol<-mapIds(org.Hs.eg.db,keys=dex_filt_ensembl,column="SYMBOL",keytype = "ENSEMBL",multiVals = "first")
dev.off()

# vennDiagram DEU, DEGs and lost m6A genes --------------------------------

#Overlap DEGs, lost m6A and DEU:
allDEGs<-diffexp_genes$gene
names(allDEGs)<-diffexp_genes$gene
allDEGs<-unique(allDEGs)

names(lost_m6a_genes_wu)<-lost_m6a_genes_wu
names(lost_m6a_genes_uk)<-lost_m6a_genes_uk
names(lost_m6a_genes_sa)<-lost_m6a_genes_sa

deu_genes<-dex_filt$gene
names(deu_genes)<-deu_genes
deu_genes<-unique(deu_genes)

deu_genes_entrez <- mapIds(org.Hs.eg.db,keys = sig_sa, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")

lostm6a_genes_all<-unique(union(union(lost_m6a_genes_wu,lost_m6a_genes_uk),lost_m6a_genes_sa))
names(lostm6a_genes_all)<-lostm6a_genes_all

ggvenn(list("DEU"=deu_genes,"DEGs"=allDEGs,"lost m6A WU"=lost_m6a_genes_wu,"lost m6A SA"=lost_m6a_genes_sa,"lost m6A UK"=lost_m6a_genes_uk))

ggvenn(list("DEU"=deu_genes,"DEGs"=allDEGs,"lost m6A"=lostm6a_genes_all))

common_diffexp_genes
names(common_diffexp_genes)<-common_diffexp_genes

ggvenn(list("DEU"=deu_genes,"common DEGs"=common_diffexp_genes,"lost m6A"=lostm6a_genes_all))

intersect(intersect(deu_genes,common_diffexp_genes),lostm6a_genes_all)

# Enrichment Analysis DEU genes -------------------------------------------

goBP_deu_genes<-enrichGO(gene=deu_genes,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F)%>% gofilter(.,level = c(7:10))
goMF_deu_genes<-enrichGO(gene=deu_genes,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goCC_deu_genes<-enrichGO(gene=deu_genes,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
kegg_deu_genes<-enrichKEGG(deu_genes_entrez,pvalueCutoff = 0.05)
reactome_deu_genes<-enrichPathway(deu_genes_entrez,pvalueCutoff=0.05)
enrichrWP_deu_genes<-enrichr(deu_genes,"WikiPathways_2019_Human")
enrichrCov2sets_deu_genes<-enrichr(deu_genes,"COVID-19_Related_Gene_Sets")

pdf("/plots/enrichmentGO_DEU_genes_all.pdf")

mydotplot(goBP_deu_genes,showCategory = 60,x="Count",title = "GO:Biological Process \n DEU genes") + theme(text = element_text(size = 6))

mydotplot(goMF_deu_genes,showCategory = 60,x="Count",title = "GO:Molecular Function \n DEU genes") + theme(text = element_text(size = 6))

mydotplot(goCC_deu_genes,showCategory = 60,x="Count",title = "GO:Cellular Component \n DEU genes") + theme(text = element_text(size = 6))

mydotplot(kegg_deu_genes,showCategory = 60,x="Count",title = "KEGG \n DEU genes") + theme(text = element_text(size = 6))

mydotplot(reactome_deu_genes,showCategory = 60,x="Count",title = "Reactome \n DEU genes")  + theme(text = element_text(size = ))

plotEnrich(enrichrWP_deu_genes$WikiPathways_2019_Human,title="Enrichr:WikiPathways \n DEU genes")

plotEnrich(enrichrCov2sets_deu_genes$`COVID-19_Related_Gene_Sets`,title="Enrichr:SARS-Cov2 sets \n DEU genes")


dev.off()

goTables_deu<-list("GO_BP"=goBP_deu_genes %>% filter(pvalue<0.05),
                   "GO_MF"=goMF_deu_genes %>% filter(pvalue<0.05),
                   "GO_CC"=goCC_deu_genes %>% filter(pvalue<0.05),
                   "KEGG"=kegg_deu_genes %>% filter(pvalue<0.05),
                   "Reactome"=reactome_deu_genes %>% filter(pvalue<0.05),
                   "Enrichr_WikiPaths"=enrichrWP_deu_genes$WikiPathways_2019_Human %>% filter(P.value<0.05),
                   "Enrichr_SARS-Cov2_sets"=enrichrCov2sets_deu_genes$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05))

openxlsx::write.xlsx(goTables_deu,file="enrichmentAnalysisTable_DEUgenes.xls")


# Intersection of DEU genes with m6A peaks in each condition: -------------


all_exonic_m6a_intersections<- fread("/reference_genomes/bedtools_intersect_Chlorocebus_sabaeus.ChlSab1.1.104.chr_ensembl_flattened_DEXseq_narrowPeaks_vero_wu_uk_sa_EXONIC_PART.tsv",header=F,sep="\t")

colnames(all_exonic_m6a_intersections) <- c(paste0(c("chr","start","end","id","score","strand","source","type","phase","features"),"_annot"),
paste0(c("condition","chr","start","end","peak_name","score","strand","signalValue","-log10pValue","-log10qValue","pointSource"),"_peaks"))


all_exonic_m6a_intersections <- all_exonic_m6a_intersections %>%
  separate(features_annot,into = c("gene_id","transcript_id","exon"),sep = ";") %>%
  mutate(gene_id=gsub("(gene_id |\")","",gene_id,perl = T),
         transcript_id=gsub("(transcripts |\")","",transcript_id,perl = T),
         exon=gsub("(exonic_part_number\\s+|\")","",exon,perl = T)) %>% 
          mutate(exon=gsub("\\s+","",exon,perl = T), gene_id=gsub("\\+.*","",gene_id), exon=paste0("E",exon))

#How many DEU genes contain with m6a peaks at the exon level? (join by gene name and exon name):
dex_filt_all_exonic_m6a_intersections <- inner_join(dex_filt,all_exonic_m6a_intersections,by=c("groupID"="gene_id","featureID"="exon"))

dex_filt_all_exonic_m6a_intersections %>% filter(!is.na(peak_name_peaks)) %>%
  dplyr::select(groupID) %>% distinct() %>% dim()

dex_filt_all_exonic_m6a_intersections %>% dplyr::select(condition_peaks,peak_name_peaks) %>% distinct() %>% group_by(condition_peaks) %>% dplyr::count()

left_join(dex_filt_all_exonic_m6a_intersections,peaks_table %>% dplyr::select(name_macs2,short_annotation),by=c("peak_name_peaks"="name_macs2")) %>% dplyr::count(short_annotation)

all_exonic_m6a_DEU_vero <- dex_filt_all_exonic_m6a_intersections %>% filter(condition_peaks=="Vero") %>% dplyr::select(gene) %>% distinct() %>% unlist()
names(all_exonic_m6a_DEU_vero) <- all_exonic_m6a_DEU_vero



#Impute empty gene names:
map_ensembl_to_entrezid<-getBM(filters ="ensembl_gene_id" , attributes = c("ensembl_gene_id","entrezgene_accession","hgnc_symbol","uniprot_gn_symbol","external_gene_name"),values =dex_filt$groupID ,mart = mart)
dex_filt<- left_join(dex_filt,map_ensembl_to_entrezid %>% dplyr::select(ensembl_gene_id,entrezgene_accession),by=c("groupID"="ensembl_gene_id")) %>% mutate(imputed_gene_name=if_else(entrezgene_accession=="",gene,entrezgene_accession)) %>% distinct()
#annotated_gene_symbols_df<- dex_filt_gene_symbol %>% mutate(gene=if_else(hgnc_symbol=="",uniprot_gn_symbol,hgnc_symbol)) %>% mutate(gene=if_else(gene=="",ensembl_gene_id,gene,ensembl_gene_id)) 

deuGenes_with_peaks_annotation<- left_join(dex_filt,peaks_table,by=c("imputed_gene_name"="peakAnn_gene_name")) %>% 
  filter(!is.na(name_macs2) &! is.na(condition) & short_annotation !="intergenic") %>% rename(condition_peak=condition)

# How many genes have m6a and DEU?
deuGenes_with_peaks_annotation %>% dplyr::select(groupID,condition_peak) %>% distinct() %>% group_by(condition_peak) %>% dplyr::count()

deuGenes_with_peaks_annotation %>% dplyr::select(groupID,condition_peak,name_macs2,short_annotation) %>% distinct() %>% group_by(condition_peak) %>% dplyr::count(short_annotation)



# Number of peaks per condition in DEU genes: ----------------------------

deuGenes_with_peaks_annotation %>%
  dplyr::select(condition_peak,name_macs2) %>% 
  distinct() %>% 
  group_by(condition_peak) %>% 
  dplyr::count() %>% ungroup()
#___________________Barplot Number of m6A peaks in DEU genes _______________________#
deuGenes_with_peaks_annotation %>%
  dplyr::select(condition_peak,groupID,name_macs2) %>% 
  distinct() %>% 
  group_by(condition_peak) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  ggplot(aes(fct_relevel(condition_peak,levels=c("vero","wu","uk","sa")),n,fill=condition_peak)) + 
  geom_bar(stat = "identity") +
  xlab("") + ylab("Number of peaks") + labs(title="Number of peaks in DEU genes",fill="") + theme_light() + theme(aspect.ratio = 1.5) 



# Peaks per gene by condition (DEU genes): --------------------------------


countsPeaksPerGeneDEU <- deuGenes_with_peaks_annotation %>%
  dplyr::select(condition_peak,groupID,name_macs2) %>% 
  distinct() %>% 
  group_by(condition_peak,groupID) %>% 
  dplyr::count()

countsPeaksPerGeneDEU

# Summary Peaks Per DEU Gene ----------------------------------------------
# Boxplot: Peaks per gene by condition (DEU genes) ------------------------


deuGenes_with_peaks_annotation %>%
  dplyr::select(condition_peak,groupID,name_macs2) %>% 
  distinct() %>% 
  group_by(condition_peak,groupID) %>% 
  dplyr::count() %>% ungroup() %>% 
  ggplot(aes(fct_relevel(condition_peak,levels=c("vero","wu","uk","sa")),n,color=condition_peak)) + 
  geom_point(alpha=0.7,position = position_jitter(height = 0.05), size=0.75) +
  geom_boxplot(alpha=0.5,notch = T, notchwidth = 0.5, varwidth = T,outlier.size = 0.5, colour="black", width=0.5)  +
  #stat_summary(fun.y = "mean", color="red") +
  theme_light() + 
  theme(aspect.ratio = 1, panel.background = element_blank(), panel.grid = element_blank(),legend.position = "none") +
  ylab("Number of peaks per gene") +
  xlab("") + 
  labs(fill="",title="Peaks per gene",subtitle ="DEU genes") +
  scale_x_discrete(labels=c("vero"="Non-infected","wu"="B.1","uk"="B.1.1.7","sa"="B.1.351")) +
  scale_y_continuous(breaks = seq(0,15), labels=seq(0,15)) + ylim(0,10) ->plotPeaksPerGene_DEU
plotPeaksPerGene_DEU

ggsave("/plots/peaksPerGene_DEUGenes.pdf",plot = plotPeaksPerGene_DEU,device = "pdf", width = 6,height = 8,dpi = 300)

  
#scale_y_continuous(breaks=c(0:10),labels=as.character(c(0:10)), limits = c(0.001,10), expand = c(0.001,0.5)) +


#Number of peaks per DEU gene:
  peaksPerDEUGene<- deuGenes_with_peaks_annotation %>%
    dplyr::select(condition_peak,groupID,name_macs2) %>% 
    distinct() %>% 
    group_by(condition_peak,groupID) %>% 
    dplyr::count() %>% ungroup()


# Summary at GENE level per condition: ------------------------------------


  numDEUGenes_with_m6a<- deuGenes_with_peaks_annotation %>%
    dplyr::select(condition_peak,groupID) %>% 
    distinct() %>% group_by(condition_peak) %>% dplyr::count() %>% as.data.frame()
  
  numDEUs_with_m6a_vero<-numDEUGenes_with_m6a$n[which(numDEUGenes_with_m6a$condition_peak=="vero")]
  numDEUs_with_m6a_wu<-numDEUGenes_with_m6a$n[which(numDEUGenes_with_m6a$condition_peak=="wu")]
  numDEUs_with_m6a_uk<-numDEUGenes_with_m6a$n[which(numDEUGenes_with_m6a$condition_peak=="uk")]
  numDEUs_with_m6a_sa<-numDEUGenes_with_m6a$n[which(numDEUGenes_with_m6a$condition_peak=="sa")]
  
  numPeaks_DEUgenes<-deuGenes_with_peaks_annotation %>%
    dplyr::select(condition_peak,groupID,name_macs2) %>% 
    distinct() %>% 
    group_by(condition_peak) %>% 
    dplyr::count()

  #summary_peaks_in_DEU<- numPeaks_DEUgenes$n / c(numDEUs_with_m6a_sa,numDEUs_with_m6a_uk,numDEUs_with_m6a_vero,numDEUs_with_m6a_wu)
  

# Table Peaks per gene summary, DEU genes ---------------------------------

  summaryPeaksPerDEUGene <- peaksPerDEUGene %>% group_by(condition_peak) %>% summarise(mean_num_peaks_per_gene=mean(n), median(n))
  summaryPeaksPerDEUGene$condition_peak <-factor(summaryPeaksPerDEUGene$condition_peak,levels = c("vero","wu","uk","sa"))
  summaryPeaksPerDEUGene %>% arrange(condition_peak) %>% kbl() %>% kable_classic(full_width=F)


# Number of Peaks by Annotation, DEU genes --------------------------------

  #Where are peaks from DEU genes located (exon, intron, etc):
  numPeaksByAnnotationDEUGenes<- deuGenes_with_peaks_annotation %>%
    dplyr::select(condition_peak,groupID,name_macs2, short_annotation) %>% 
    distinct() %>% 
    group_by(condition_peak, short_annotation) %>% 
    dplyr::count() %>% ungroup()
  
  numPeaksByAnnotationDEUGenes %>% 
    ggplot(aes(fct_relevel(condition_peak,levels=c("vero","wu","uk","sa")),n,fill=short_annotation)) +
    geom_bar(stat = "identity", position = position_dodge()) + 
    xlab("Condition") + ylab("Number of peaks per annotaion") + labs(title = "Number of peaks by annotation", subtitle = "DEU genes")
  
  #Number of DEU genes with direct intersection with m6A peaks:
  dex_filt_all_exonic_m6a_intersections %>% dplyr::select(condition_peaks,groupID) %>% group_by(condition_peaks) %>% dplyr::count()
  
  #Number of peaks in DEU genes with direct intersection of m6A:
  dex_filt_all_exonic_m6a_intersections %>% dplyr::select(condition_peaks,groupID,peak_name_peaks) %>% group_by(condition_peaks,groupID) %>% dplyr::count()
  
#DEU genes having m6A modification directly at the predicted DEU exons
#Number of peaks per gene:
  peaksPerDEUIntersectingExonm6A<-dex_filt_all_exonic_m6a_intersections %>% 
    dplyr::select(condition_peaks,groupID,peak_name_peaks) %>%
    distinct() %>%
    group_by(condition_peaks,groupID) %>% 
    dplyr::count() %>% ungroup()
  
  
# Table PeaksPerGene, DEUs with intersecting m6A at exon ------------------

  peaksPerDEUIntersectingExonm6A$condition_peaks <- factor(peaksPerDEUIntersectingExonm6A$condition_peaks,levels = c("Vero","WU","UK","SA"))
  peaksPerDEUIntersectingExonm6A %>% group_by(condition_peaks) %>% summarise(mean_num_peaks_per_gene=mean(n), median(n))
  
#Number of peaks per exon:
peaksPerExonDEUIntersectingExonm6A<-dex_filt_all_exonic_m6a_intersections %>% 
  dplyr::select(condition_peaks,groupID,featureID,peak_name_peaks) %>%
    distinct() %>%
    group_by(condition_peaks,groupID,featureID) %>% 
    dplyr::count() %>% ungroup()

# Table Peaks Per Exon m6A intersecting Exon -----------------------------------

  peaksPerExonDEUIntersectingExonm6A$condition_peaks <- factor(peaksPerExonDEUIntersectingExonm6A$condition_peaks,levels = c("Vero","WU","UK","SA"))
  peaksPerExonDEUIntersectingExonm6A %>% group_by(condition_peaks) %>% 
    summarise(mean_num_peaks_per_gene=mean(n), median(n)) %>%
    arrange(condition_peaks) %>% kbl() %>% kable_classic(full_width=F)
  

  

# Boxplot peaks per gene in DEU genes intersecting m6A at exons -----------

  peaksPerGene_intersectingm6AtExon<- dex_filt_all_exonic_m6a_intersections %>%
  dplyr::select(condition_peaks,groupID,peak_name_peaks) %>% 
  distinct() %>% 
  group_by(condition_peaks,groupID) %>% 
  dplyr::count() %>%
  ungroup() 
  
  # Table Peaks Per Gene, DEU genes with Intersecting m6A at Exon -----------
  
  peaksPerGene_intersectingm6AtExon <- peaksPerGene_intersectingm6AtExon %>% group_by(condition_peaks) %>% summarise(mean_num_peaks_per_gene=mean(n), median(n))
  peaksPerGene_intersectingm6AtExon$condition_peaks <-factor(summaryPeaksPerNonDEUGene$condition_peak,levels = c("Vero","WU","UK","SA"))
  peaksPerGene_intersectingm6AtExon %>% arrange(condition_peaks) %>% kbl() %>% kable_classic(full_width=F)
  
  
  peaksPerGene_intersectingm6AtExon %>% 
  ggplot(aes(fct_relevel(condition_peaks,levels=c("Vero","WU","UK","SA")),n,color=condition_peaks)) + 
  geom_point(alpha=0.1,position = position_jitterdodge()) +
  geom_boxplot(alpha=0.5,notch = F, notchwidth = 0.75, varwidth = T,outlier.size = 0.5, colour="black")  +
  stat_summary(fun = "mean", color="red")+
  theme_light() + theme(aspect.ratio = 1.5, panel.background = element_blank(), panel.grid = element_blank()) +
  ylab("Number of peaks per gene") + xlab("") + labs(fill="", title="Peaks per gene",subtitle = "m6A overlapping with DEU exons") -> plot_numPeaksGeneDEU_intersecting_m6A
  plot_numPeaksGeneDEU_intersecting_m6A

  # Boxplot peaks per Exon in DEU genes intersecting m6A at exons -----------
  
  dex_filt_all_exonic_m6a_intersections %>%
    dplyr::select(condition_peaks,groupID,featureID,peak_name_peaks) %>% 
    distinct() %>% 
    group_by(condition_peaks,groupID,featureID) %>% 
    dplyr::count() %>%
    ungroup() %>% 
    ggplot(aes(fct_relevel(condition_peaks,levels=c("vero","WU","UK","SA")),n,color=condition_peaks)) + 
    geom_point(alpha=0.1,position = position_jitterdodge()) +
    geom_boxplot(alpha=0.5,notch = F, notchwidth = 0.75, varwidth = T,outlier.size = 0.5, colour="black")  +
    stat_summary(fun.y = "mean", color="red") +
    theme_light() + theme(aspect.ratio = 1.5, panel.background = element_blank(), panel.grid = element_blank()) +
    ylab("Number of peaks per exon") + xlab("") + labs(fill="", title="Peaks per exon",subtitle = "m6A overlapping with DEU exons")





# Peaks Per Gene , Non-DEU genes ------------------------------------------


peaksPerGene_nonDEU<- allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
  distinct() %>% 
  filter(!(deg_gene %in% unique(dex_filt$imputed_gene_name))) %>% 
  group_by(condition_peak,deg_gene) %>% 
  dplyr::count() %>% 
  ungroup() 

# Table, summary peaks per gene, Non-DEU genes ----------------------------

summaryPeaksPerNonDEUGene <- peaksPerGene_nonDEU %>% group_by(condition_peak) %>% summarise(mean_num_peaks_per_gene=mean(n), median(n))
summaryPeaksPerNonDEUGene$condition_peak <-factor(summaryPeaksPerNonDEUGene$condition_peak,levels = c("vero","WU","UK","SA"))
summaryPeaksPerNonDEUGene %>% arrange(condition_peak) %>% kbl() %>% kable_classic(full_width=F)


# Boxplot Peaks Per Gene, Non-DEU genes -----------------------------------

peaksPerGene_nonDEU %>% 
  ggplot(aes(fct_relevel(condition_peak,levels=c("vero","WU","UK","SA")),n,color=condition_peak)) + 
  #stat_boxplot(geom = "errorbar",width=0.2, size=0.5, colour="black") +
  geom_point(alpha=0.1,position = position_jitter(), size=0.1) +
  geom_boxplot(alpha=0.7,notch = T, notchwidth = 0.5, varwidth = T,outlier.size = 0.1,outlier.alpha = 0.5, colour="black", width=0.5)  +
  #stat_summary(fun.y = "mean", color="red") +
  theme_light() + theme(aspect.ratio = 1.5, panel.background = element_blank(), panel.grid = element_blank(), legend.position = "none") +
  ylab("Number of peaks per gene") + xlab("") + labs(color="",fill="", title = "Number of peaks per gene",subtitle = "Non-DEU genes") +
  scale_x_discrete(labels=c("vero"="Non-infected","WU"="B.1","UK"="B.1.1.7","SA"="B.1.351")) +
  scale_y_discrete(limits=c(0,15), labels=seq(0,15)) + ylim(0,15)->plotPeaksPerGene_nonDEU
plotPeaksPerGene_nonDEU

ggsave("/plots/peaksPerGene_nonDEUGenes.pdf",plot = plotPeaksPerGene_nonDEU,device = "pdf", width = 6,height = 8,dpi = 300)


#Counting peaks in randomly selected genes from the total with a group size equal to the number of DEU genes:
numberDEUgenes_vero<- peaksPerDEUGene %>% filter(condition_peak=="vero") %>% dplyr::select(groupID) %>% distinct() %>% unlist() %>% length()
numberDEUgenes_wu<- peaksPerDEUGene %>% filter(condition_peak=="wu") %>% dplyr::select(groupID) %>% distinct() %>% unlist() %>% length()
numberDEUgenes_uk<- peaksPerDEUGene %>% filter(condition_peak=="uk") %>% dplyr::select(groupID) %>% distinct() %>% unlist() %>% length()
numberDEUgenes_sa<- peaksPerDEUGene %>% filter(condition_peak=="sa") %>% dplyr::select(groupID) %>% distinct() %>% unlist() %>% length()

randomGenes_vero<-sample( allgenes_and_peaks_retained_lost_gained %>% 
                            filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="vero")  %>% 
                      dplyr::select(deg_gene) %>% 
                      distinct() %>% unlist(),size=numberDEUgenes_vero, replace = F)

randomGenes_wu<-sample( allgenes_and_peaks_retained_lost_gained %>% 
                            filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="WU")  %>% 
                            dplyr::select(deg_gene) %>% 
                            distinct() %>% unlist(),size=numberDEUgenes_wu, replace = F)
randomGenes_uk<-sample( allgenes_and_peaks_retained_lost_gained %>% 
                            filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="UK")  %>% 
                            dplyr::select(deg_gene) %>% 
                            distinct() %>% unlist(),size=numberDEUgenes_uk, replace = F)
randomGenes_sa<-sample( allgenes_and_peaks_retained_lost_gained %>% 
                            filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="SA")  %>% 
                            dplyr::select(deg_gene) %>% 
                            distinct() %>% unlist(),size=numberDEUgenes_sa, replace = F)


# Peaks  per Gene, Random Genes -------------------------------------------


randomPeakCounts_vero <- allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="vero")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
  distinct() %>% 
  filter(deg_gene %in% unique(randomGenes_vero)) %>% 
  group_by(condition_peak,deg_gene) %>% 
  dplyr::count() %>% ungroup()

randomPeakCounts_wu<-allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="WU")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
  distinct() %>% 
  filter(deg_gene %in% unique(randomGenes_wu)) %>% 
  group_by(condition_peak,deg_gene) %>% 
  dplyr::count() %>% ungroup()

randomPeakCounts_uk<-allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="UK")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
  distinct() %>% 
  filter(deg_gene %in% unique(randomGenes_uk)) %>% 
  group_by(condition_peak,deg_gene) %>% 
  dplyr::count() %>% ungroup()

randomPeakCounts_sa<-allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="SA")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
  distinct() %>% 
  filter(deg_gene %in% unique(randomGenes_sa)) %>% 
  group_by(condition_peak,deg_gene) %>% 
  dplyr::count() %>% ungroup()

randomPeakCounts <-rbind(randomPeakCounts_vero,randomPeakCounts_wu,randomPeakCounts_uk,randomPeakCounts_sa)

# Table, summary peaks per gene, Random Genes -----------------------------

summaryrandomPeakCounts <- randomPeakCounts %>% group_by(condition_peak) %>% summarise(mean_num_peaks_per_gene=mean(n), median(n))
summaryrandomPeakCounts$condition_peak <-factor(summaryrandomPeakCounts$condition_peak,levels = c("vero","WU","UK","SA"))
summaryrandomPeakCounts %>% arrange(condition_peak) %>% kbl() %>% kable_classic(full_width=F)


# Boxplot Peaks Per Gene, Random Genes ------------------------------------


randomPeakCounts %>% 
  ggplot(aes(fct_relevel(condition_peak,levels=c("vero","WU","UK","SA")),n,color=condition_peak)) + 
  #stat_boxplot(geom = "errorbar",width=0.2, size=0.5) +
  geom_point(alpha=0.7,position = position_jitter(height = 0.05), size=0.75) +
  geom_boxplot(notch = T, notchwidth = 0.5, varwidth = T,outlier.size = 0.5, colour="black", width=0.5, alpha=0.7)  +
  #stat_summary(fun.y = "mean", color="red") +
  theme_light() + theme(aspect.ratio = 1.5, panel.background = element_blank(), panel.grid = element_blank(), legend.position = "none") +
  ylab("Number of peaks per gene") + xlab("") + labs(fill="", title = "Number of peaks per gene",subtitle = "Random genes") +
  scale_x_discrete(labels=c("vero"="Non-infected","WU"="B.1","UK"="B.1.1.7","SA"="B.1.351")) +
  scale_y_continuous(breaks=seq(0,15), labels=seq(0,15)) + ylim(0,15) -> plotPeaksPerGene_random
plotPeaksPerGene_random


ggsave("/plots/peaksPerGene_random.pdf",plot = plotPeaksPerGene_random,device = "pdf", width = 6,height = 8,dpi = 300)


# Peaks per gene, DEU vs Random -------------------------------------------
randomSampCountsDEU <- randomPeakCounts %>% filter(condition_peak=="vero") %>% mutate(type="Random")
deu_and_random<- rbind(countsPeaksPerGeneDEU %>% ungroup() %>% rename(deg_gene=groupID) %>% filter(condition_peak=="vero") %>% mutate(type="DEU"),randomSampCountsDEU)

deu_and_random %>% ggplot(aes(fct_relevel(type,levels=c("vero","wu","uk","sa")),n)) + 
  geom_point(alpha=0.7,position = position_jitter(height = 0.05), size=0.15, color="steelblue") +
  geom_boxplot(alpha=0.5,notch = T, notchwidth = 0.5, varwidth = T,outlier.size = 0.5, colour="black", width=0.5)  +
  #stat_summary(fun.y = "mean", color="red") +
  theme_light() + 
  theme(aspect.ratio = 1.25, panel.background = element_blank(), panel.grid = element_blank(),legend.position = "none") +
  ylab("Number of peaks per gene") +
  xlab("") + 
  labs(fill="",title="Peaks per gene",subtitle ="DEU genes") +
  scale_x_discrete(labels=c("vero"="Non-infected","wu"="B.1","uk"="B.1.1.7","sa"="B.1.351")) +
  scale_y_continuous(breaks = seq(0,15), labels=seq(0,15)) + ylim(0,10) ->plotPeaksPerGene_DEU_random

plotPeaksPerGene_DEU_random

ggsave("/plots/peaksPerGene_DEU_nonInfected_vs_random_AND_peaksPerGeneDEU_conditions.pdf",plot_grid(plotPeaksPerGene_DEU_random,plotPeaksPerGene_DEU),width = 10,height = 8,dpi = 300)


# Null distribution -------------------------------------------------------

peaksPerGeneVero <- allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="vero")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
  distinct() %>% group_by(condition_peak,deg_gene) %>% 
  dplyr::count() %>% ungroup()

point_estimate_mean <- countsPeaksPerGeneDEU %>% filter(condition_peak=="vero") %>% 
  ungroup() %>% specify(response = n) %>% calculate(stat="mean")

point_estimate_median <- countsPeaksPerGeneDEU %>% filter(condition_peak=="vero") %>% 
  ungroup() %>% specify(response = n) %>% calculate(stat="median")

null_dist <- peaksPerGeneVero %>% filter(condition_peak=="vero") %>% 
  ungroup() %>% specify(response = n) %>% 
  hypothesise(null = "point", med = 2) %>%
  generate(reps = 1000, type = "bootstrap") 

#Run once:
  getRandomGenes <- function(df,n, size){
    set.seed(123)
    sampdf <- NULL
    
    for(i in 1:n){
      
      randomGenes<-sample( df %>% filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="vero")  %>% 
                                  dplyr::select(deg_gene) %>% 
                                  distinct() %>% unlist(),size=size, replace = F)
      
      df1 <- df %>% 
        filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic" & condition_peak=="vero")  %>% 
        dplyr::select(condition_peak,deg_gene,name_macs2) %>% 
        distinct() %>% 
        filter(deg_gene %in% unique(randomGenes)) %>% 
        group_by(condition_peak,deg_gene) %>% 
        dplyr::count() %>% ungroup() %>% 
        mutate(replicate=i)
      
      sampdf <-rbind(sampdf,df1)
      df1<-NULL
    }
    
    return(sampdf)
  }

  
  #randomPeakCounts_vero_samp1000 <- getRandomGenes(df=allgenes_and_peaks_retained_lost_gained,n=1000,size=numberDEUgenes_vero)
  
  #save(randomPeakCounts_vero_samp1000,file="/plots/tables/randomSampling1000rep_DEUgenestest_vero.Rda")
  #load(file="/plots/tables/randomSampling1000rep_DEUgenestest_vero.Rda") 
  
  
  #mean_random_dist <- randomPeakCounts_vero_samp1000 %>% group_by(replicate) %>% summarise(mean=mean(n)) %>% summarise(mean(mean)) %>% unlist()
  
  pdf("/plots/numPeaksDEUgenes_vs_random.pdf", width = 4, height =3)

  randomPeakCounts_vero_samp1000 %>% group_by(replicate) %>% summarise(mean=mean(n)) %>% 
    ggplot(aes(mean)) + geom_histogram(bins = 100) + 
    geom_vline(xintercept = unlist(point_estimate_mean), color="red") +
    geom_vline(xintercept = mean_random_dist, color="blue") +
    theme_bw(base_size = 11) + 
    xlab("Mean number of peaks") +
    ylab("Count") +
    annotate(geom = "text", x=2.2,y=75,label="t-test \n p < 2.2e-16") +
    scale_y_continuous(breaks=seq(0,80,by=10)) +
    scale_x_continuous(breaks=seq(1.5,3,by=0.25))
  
  dev.off()
  
  randomPeakCounts_vero_samp1000 %>% group_by(replicate) %>% summarise(median=median(n)) %>% 
    ggplot(aes(median)) + geom_histogram(bins = 100) + geom_vline(xintercept = unlist(point_estimate_median), color="red") +
    theme_bw(base_size = 11) + xlab("Median number of peaks")
  
  t.test(randomPeakCounts_vero_samp1000 %>% group_by(replicate) %>% summarise(mean=mean(n)) %>% dplyr::select(mean) %>% unlist(), mu=point_estimate_mean$stat)


#Permutation test
#############
permutation.test <- function(treatment, outcome, n){
  distribution=c()
  result=0
  for(i in 1:n){
    distribution[i]=diff(by(outcome, sample(treatment, length(treatment), FALSE), mean))
  }
  result=sum(abs(distribution) >= abs(original+1))/(n+1)
  return(list(result, distribution))
}

test1 <- permutation.test($n, randomPeakCounts$n, 10000)
hist(test1[[2]], breaks=50, col='grey', main="Permutation Distribution", las=1, xlab='')
abline(v=original, lwd=3, col="red")

test1[[1]]




################
#m6A fold enrichment in Non-DEU genes
allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2, pileup_macs2,  fold_enrichment_macs2) %>% 
  distinct() %>% 
  filter(!(deg_gene %in% unique(dex_filt$gene))) %>% 
  group_by(condition_peak) %>% 
  ggplot(aes(fct_relevel(condition_peak,levels=c("vero","WU","UK","SA")),log2(fold_enrichment_macs2),color=condition_peak)) + 
  geom_violin(alpha=0.25) + geom_boxplot(alpha=0.5)  +
  ylab("log2 Fold m6A Enrichment") + xlab("") + labs(fill="", title = "m6A Fold enrichment",subtitle = "Non-DEU genes") +
  theme_light() + theme(aspect.ratio = 1.5)


allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2, pileup_macs2, fold_enrichment_macs2) %>% 
  distinct() %>% filter(deg_gene %in% unique(dex_filt$imputed_gene_name)) %>% 
  group_by(condition_peak) %>% 
  ggplot(aes(fct_relevel(condition_peak,levels=c("vero","WU","UK","SA")), log2(fold_enrichment_macs2),color=condition_peak)) + 
  geom_violin(alpha=0.25) + geom_boxplot(alpha=0.5)  +
  ylab("log2 Fold m6A Enrichment") + xlab("") + labs(fill="", title = "m6A Fold enrichment",subtitle = "DEU genes") +
  theme_light() + theme(aspect.ratio = 1.5)



allgenes_and_peaks_retained_lost_gained %>% 
  filter(!is.na(condition_peak) &!is.na(name_macs2) & short_annotation !="intergenic")  %>% 
  dplyr::select(condition_peak,deg_gene,name_macs2, pileup_macs2, fold_enrichment_macs2) %>% 
  distinct() %>% filter(deg_gene %in% unique(deuGenes_with_peaks_annotation$imputed_gene_name)) %>% 
  group_by(condition_peak) %>% 
  ggplot(aes(fct_relevel(condition_peak,levels=c("vero","WU","UK","SA")),log2(fold_enrichment_macs2),color=condition_peak)) + 
  geom_violin(alpha=0.25) + geom_boxplot(alpha=0.5)  +
  ylab("log2 Fold m6A Enrichment") + xlab("") + labs(fill="", title = "m6A Fold Enrichment",subtitle = "DEU genes overlapping retained m6A peaks") +
  theme_light() + theme(aspect.ratio = 1.5)


# Selection of genes for primer design:
##############

# 1. Criteria: DEU genes pertaining to common DEGs which have m6A modification directly in the DEU predicted exons:

top10_deu_commonDEG_intersecting_exon_m6a <- dex_filt %>% 
  filter((gene %in% common_diffexp_genes) & (gene %in%  dex_filt_all_exonic_m6a_intersections$gene)) %>%
  dplyr::select(-transcripts) %>% 
  pivot_longer(names_to = "contrast",values_to = "log2FC",c("log2fold_WU_Vero","log2fold_UK_Vero","log2fold_SA_Vero")) %>%
  group_by(contrast) %>% 
  arrange(dplyr::desc(abs(log2FC))) %>% 
  top_n(10) %>% 
  arrange(contrast) %>%
  as.data.frame()

openxlsx::write.xlsx(list(top10_deu_commonDEG_intersecting_exon_m6a),"/plots/manuscript_figures/tables/top_deu_commonDEG_intersecting_exon_m6a.xls")




# Analyze shared DEU genes across strains:
# How many DEUs (log2FC > 1) by strain?:
dex_filt_long <- dex_filt %>% 
  pivot_longer(names_to = "contrast",values_to = "log2FC",c("log2fold_WU_Vero","log2fold_UK_Vero","log2fold_SA_Vero")) 

#What is the mean and median values of log2FC across strains?:
dex_filt_long %>% summarise(mean(log2FC),median(log2FC))

threshold<-log2(1.1)

dex_filt_long %>%
  group_by(contrast) %>%
  filter(abs(log2FC) > threshold) %>% 
  dplyr::select(contrast,gene) %>% 
  distinct() %>% 
  dplyr::count() %>% 
  kbl() %>%
  kable_classic(full_width = F, html_font = "Cambria")



deu_genes_wu_abslogFC <- dex_filt %>% 
  pivot_longer(names_to = "contrast",values_to = "log2FC",c("log2fold_WU_Vero","log2fold_UK_Vero","log2fold_SA_Vero")) %>%
  filter(contrast=="log2fold_WU_Vero") %>% 
  filter(abs(log2FC) > threshold) %>% 
  dplyr::select(gene) %>% 
  distinct() %>% 
  unlist()

deu_genes_uk_abslogFC <- dex_filt %>% 
  pivot_longer(names_to = "contrast",values_to = "log2FC",c("log2fold_WU_Vero","log2fold_UK_Vero","log2fold_SA_Vero")) %>%
  filter(contrast=="log2fold_UK_Vero") %>% 
  filter(abs(log2FC) > threshold) %>% 
  dplyr::select(gene) %>% 
  distinct() %>% 
  unlist()

deu_genes_sa_abslogFC <- dex_filt %>% 
  pivot_longer(names_to = "contrast",values_to = "log2FC",c("log2fold_WU_Vero","log2fold_UK_Vero","log2fold_SA_Vero")) %>%
  filter(contrast=="log2fold_SA_Vero") %>% 
  filter(abs(log2FC) > threshold) %>% 
  dplyr::select(gene) %>% 
  distinct() %>% 
  unlist()
  
commonDEUs_above_thres<-intersect(intersect(deu_genes_wu_abslogFC,deu_genes_uk_abslogFC),deu_genes_sa_abslogFC)

DEUs_above_thres<-union(union(deu_genes_wu_abslogFC,deu_genes_uk_abslogFC),deu_genes_sa_abslogFC)
pdf("/plots/manuscript_figures/vennDiagram_DEUs_across_strains_threshol10percent.pdf")

ggvenn(list("DEU WU"=deu_genes_wu_abslogFC,"DEU UK"=deu_genes_uk_abslogFC,"DEU SA" =deu_genes_sa_abslogFC)) + labs(title="DEU genes log2FC above threshold (10% fold change)")

dev.off()
openxlsx::write.xlsx(list("Common DEUs"=commonDEUs_above_thres,"DEU WU"=deu_genes_wu_abslogFC,"DEU UK"=deu_genes_uk_abslogFC,"DEU SA" =deu_genes_sa_abslogFC),file = "/plots/manuscript_figures/tables/DEU_genes_list_common_across_strains.xls")



##### Save gene coordinates for scanning in RNA binding databases:
deuGenes_bed <-flatGTF %>% 
  filter(IDtype_gtf=="aggregate_gene" & features_gtf %in% dex_filt$groupID) %>%
  dplyr::select(chr_gtf,start_gtf,end_gtf,features_gtf,score_gtf,strand_gtf) %>%
  as.data.frame()

fwrite(deuGenes_bed,file = "/plots/manuscript_figures/tables/DEUgenes_coordinates.bed", col.names = F,sep = "\t")


#geneList<-goTables_DEGs_all  %>% arrange(dplyr::desc(GeneRatio)) %>% filter(grepl("(nuclear|mRNA|catabolic|SRP|virus|viral|transcr|metabolic)",Description,ignore.case = T,perl=T)) %>% dplyr::select(geneID) %>% separate_rows(geneID,sep="/") %>% distinct() %>% unlist()

#geneList<-goTables_DEGs_all  %>% arrange(dplyr::desc(GeneRatio)) %>% filter(grepl("(nuclear|mRNA catabolic | RNA |SRP|virus|viral|transcr|metabolic| cilium|histone)",Description,ignore.case = T,perl=T)) %>% dplyr::select(Description,geneID, GeneRatio) %>% separate_rows(geneID,sep="/") %>% distinct()

geneList_wu<-go_bp_wu@result %>% arrange(dplyr::desc(GeneRatio)) %>% filter(grepl("(nuclear|mRNA catabolic | RNA |SRP|virus|viral|transcr|metabolic| cilium|histone)",Description,ignore.case = T,perl=T)) %>% dplyr::select(Description,geneID, GeneRatio) %>% separate_rows(geneID,sep="/") %>% distinct()
geneList_uk<-go_bp_uk@result %>% arrange(dplyr::desc(GeneRatio)) %>% filter(grepl("(nuclear|mRNA catabolic | RNA |SRP|virus|viral|transcr|metabolic| cilium|histone)",Description,ignore.case = T,perl=T)) %>% dplyr::select(Description,geneID, GeneRatio) %>% separate_rows(geneID,sep="/") %>% distinct()
geneList_sa<-go_bp_sa@result %>% arrange(dplyr::desc(GeneRatio)) %>% filter(grepl("(nuclear|mRNA catabolic | RNA |SRP|virus|viral|transcr|metabolic| cilium|histone)",Description,ignore.case = T,perl=T)) %>% dplyr::select(Description,geneID, GeneRatio) %>% separate_rows(geneID,sep="/") %>% distinct()

geneList_merge<-inner_join(allgenes_and_peaks_retained_lost_gained, geneList,by=c("deg_gene"="geneID")) 

geneList_merge %>% filter(deg_gene %in% lost_m6a_genes_wu & deg_contrast_numerator=="wu" & deg_gene %in% lost_m6a_genes_wu & diff_expressed==T & abs(deg_log2FoldChange) > 1 & deg_padj<0.01)


geneList_merge %>% filter(deg_gene %in% lost_m6a_genes_wu & deg_contrast_numerator=="wu" & deg_gene %in% lost_m6a_genes_wu & diff_expressed==T & abs(deg_log2FoldChange) > 1 & deg_padj<0.01)


lostPeaksTable<- allgenes_and_peaks_retained_lost_gained %>% filter(name_macs2 %in% lost_after_wu$name & short_annotation != "intergenic" & condition_peak=="vero") %>% dplyr::select(deg_gene,name_macs2,pileup_macs2,chr_macs2,start_macs2,end_macs2,fold_enrichment_macs2,condition_peak,short_annotation) %>% distinct()


inner_join(res_all_lfcShrink,peaks_table %>% filter(name_macs2 %in% lost_after_wu$name),by=c("gene"="peakAnn_gene_name")) %>%
  filter(!grepl("intergenic",short_annotation) & contrast=="wu_vs_vero" & abs(log2FoldChange)>1 & padj<0.01 &! is.na(padj)) %>% dplyr::select(gene,log2FoldChange,padj,contrast,chr_macs2,start_macs2,end_macs2,abs_summit_macs2,pileup_macs2,fold_enrichment_macs2,`-log10(qvalue)_macs2`,name_macs2,condition,orientation,short_annotation) 



infoLostAfter_wu <-inner_join(res_all_lfcShrink,peaks_table %>% filter(name_macs2 %in% lost_after_wu$name ),by=c("gene"="peakAnn_gene_name")) %>%
  filter(!grepl("(intergenic| intron)",short_annotation,perl = T) & gene %in% lost_m6a_genes_wu & contrast=="wu_vs_vero" & abs(log2FoldChange)>1 & padj<0.01 &! is.na(padj)) %>% 
  dplyr::select(gene,log2FoldChange,padj,contrast,chr_macs2,start_macs2,end_macs2,abs_summit_macs2,pileup_macs2,fold_enrichment_macs2,`-log10(qvalue)_macs2`,name_macs2,condition,orientation,short_annotation)

infoLostAfter_uk <-inner_join(res_all_lfcShrink,peaks_table %>% filter(name_macs2 %in% lost_after_uk$name ),by=c("gene"="peakAnn_gene_name")) %>%
  filter(!grepl("intergenic",short_annotation) & gene %in% lost_m6a_genes_uk & contrast=="uk_vs_vero" & abs(log2FoldChange)>1 & padj<0.01 &! is.na(padj)) %>% 
  dplyr::select(gene,log2FoldChange,padj,contrast,chr_macs2,start_macs2,end_macs2,abs_summit_macs2,pileup_macs2,fold_enrichment_macs2,`-log10(qvalue)_macs2`,name_macs2,condition,orientation,short_annotation)

infoLostAfter_sa <-inner_join(res_all_lfcShrink,peaks_table %>% filter(name_macs2 %in% lost_after_sa$name ),by=c("gene"="peakAnn_gene_name")) %>%
  filter(!grepl("intergenic",short_annotation) & gene %in% lost_m6a_genes_uk & contrast=="sa_vs_vero" & abs(log2FoldChange)>1 & padj<0.01 &! is.na(padj)) %>% 
  dplyr::select(gene,log2FoldChange,padj,contrast,chr_macs2,start_macs2,end_macs2,abs_summit_macs2,pileup_macs2,fold_enrichment_macs2,`-log10(qvalue)_macs2`,name_macs2,condition,orientation,short_annotation)


infoGainedAfter_wu <-inner_join(res_all_lfcShrink,peaks_table %>% filter(name_macs2 %in% gained_after_wu$name),by=c("gene"="peakAnn_gene_name")) %>%
  filter(!grepl("intergenic",short_annotation) & gene %in% gained_m6a_genes_wu & contrast=="wu_vs_vero" & abs(log2FoldChange)>1 & padj<0.01 &! is.na(padj)) %>% 
  dplyr::select(gene,log2FoldChange,padj,contrast,chr_macs2,start_macs2,end_macs2,abs_summit_macs2,pileup_macs2,fold_enrichment_macs2,`-log10(qvalue)_macs2`,name_macs2,condition,orientation,short_annotation)

infoGainedAfter_uk <-inner_join(res_all_lfcShrink,peaks_table %>% filter(name_macs2 %in% gained_after_uk$name),by=c("gene"="peakAnn_gene_name")) %>%
  filter(!grepl("intergenic",short_annotation) & gene %in% gained_m6a_genes_wu & contrast=="uk_vs_vero" & abs(log2FoldChange)>1 & padj<0.01 &! is.na(padj)) %>% 
  dplyr::select(gene,log2FoldChange,padj,contrast,chr_macs2,start_macs2,end_macs2,abs_summit_macs2,pileup_macs2,fold_enrichment_macs2,`-log10(qvalue)_macs2`,name_macs2,condition,orientation,short_annotation)

infoGainedAfter_sa <-inner_join(res_all_lfcShrink,peaks_table %>% filter(name_macs2 %in% gained_after_sa$name),by=c("gene"="peakAnn_gene_name")) %>%
  filter(!grepl("intergenic",short_annotation) & gene %in% gained_m6a_genes_wu & contrast=="sa_vs_vero" & abs(log2FoldChange)>1 & padj<0.01 &! is.na(padj)) %>% 
  dplyr::select(gene,log2FoldChange,padj,contrast,chr_macs2,start_macs2,end_macs2,abs_summit_macs2,pileup_macs2,fold_enrichment_macs2,`-log10(qvalue)_macs2`,name_macs2,condition,orientation,short_annotation)


infoLostAfter_wu_goTerms <- inner_join(infoLostAfter_wu,geneList_wu, by=c("gene"="geneID"))
infoLostAfter_uk_goTerms <- inner_join(infoLostAfter_uk,geneList_uk, by=c("gene"="geneID"))
infoLostAfter_sa_goTerms <- inner_join(infoLostAfter_sa,geneList_sa, by=c("gene"="geneID"))

infoGainedAfter_wu_goTerms <- inner_join(infoGainedAfter_wu,geneList_wu, by=c("gene"="geneID"))
infoGainedAfter_uk_goTerms <- inner_join(infoGainedAfter_uk,geneList_uk, by=c("gene"="geneID"))
infoGainedAfter_sa_goTerms <- inner_join(infoGainedAfter_sa,geneList_sa, by=c("gene"="geneID"))


log2FC_deu<-dex_filt %>% dplyr::select(log2fold_WU_Vero,log2fold_UK_Vero,log2fold_SA_Vero) %>% distinct()

rownames(log2FC_deu) <-dex_filt %>% 
  dplyr::select(groupID,featureID,log2fold_WU_Vero,log2fold_UK_Vero,log2fold_SA_Vero) %>%
  distinct() %>%
  mutate(label=paste0(groupID,"_",featureID)) %>%
  dplyr::select(label) %>% unlist()

infoLostAfter_wu_goTerms %>% 
  filter(grepl("splicing",Description,ignore.case = T,perl=T) &! grepl("intron",short_annotation)) %>% 
  mutate(coords=paste0(chr_macs2,":",start_macs2,"-",end_macs2)) %>% 
  arrange(dplyr::desc(log2FoldChange,fold_enrichment_macs2)) %>% 
  dplyr::select(gene,log2FoldChange,fold_enrichment_macs2,pileup_macs2,short_annotation,contrast,condition,Description,GeneRatio,orientation,coords)

#primer design
plotDEXSeq(dex_results,"ENSCSAG00000016266",legend = TRUE,cex.axis=1.2,cex=1.3,lwd=2, displayTranscripts=TRUE, norCounts=TRUE, splicing=TRUE,expression = FALSE)
flatGTF %>% filter(grepl("ENSCSAG00000012928",features_gtf)) %>% mutate(coords=paste0("chr",chr_gtf,":",start_gtf,"-",end_gtf))




# Plot: SuppFig4 example of DEU gene --------------------------------------
pdf("/plots/SuppFig4B_DEU_DEXseq_example.pdf",width = 8, height = 6)

plotDEXSeq(dex_results,"ENSCSAG00000012928",legend = TRUE,cex.axis=1.2,cex=1.3,lwd=2, displayTranscripts=F, norCounts=F, splicing=TRUE,expression = F)

dev.off()


#Correlation plots between conditions:
ggscatter(data = as.data.frame(dex_filt),x="log2fold_WU_Vero",y="log2fold_UK_Vero", add="reg.line", add.params = list(color="blue",fill="lightgray"), conf.int = T) + stat_cor(method = "spearman") + ylim(-5,5) + xlim(-5,5)

ggscatter(data = as.data.frame(dex_filt),x="log2fold_WU_Vero",y="log2fold_SA_Vero", add="reg.line", add.params = list(color="blue",fill="lightgray"), conf.int = T) + stat_cor(method = "spearman")  + ylim(-5,5) + xlim(-5,5)

ggscatter(data = as.data.frame(dex_filt),x="log2fold_SA_Vero",y="log2fold_UK_Vero", add="reg.line", add.params = list(color="blue",fill="lightgray"), conf.int = T) + stat_cor(method = "spearman")  + ylim(-5,5) + xlim(-5,5)

# Isoform switchng analyzer:
library(IsoformSwitchAnalyzeR)

flatGTF<-fread(flattenedFile,header = F,sep="\t")
colnames(flatGTF)<-paste0(c("chr","name","IDtype","start","end","score","strand","source","features"),"_gtf")

flatGTF <- flatGTF %>% filter(IDtype_gtf=="aggregate_gene") %>%
  mutate(chr_gtf=paste0("chr",chr_gtf),features_gtf=gsub("(gene_id |\"|\\+.*)","",features_gtf,perl = T)) 
#%>%  dplyr::select(chr_gtf,start_gtf,end_gtf,strand_gtf,features_gtf)

deuGenes_with_peaks_annotation<-left_join(deuGenes_with_peaks_annotation,flatGTF,by=c("groupID"="features_gtf"))

#Density plot
deuGenes_with_peaks_annotation %>%group_by(imputed_gene_name) %>% mutate(bin=(abs_summit_macs2-start_gtf)/(end_gtf-start_gtf)) %>% as.data.frame() %>% filter(bin > 0) %>% ggplot(.,aes(bin,color=condition_peak)) + geom_density() + xlim(0,1)

# extract overlapping peaks in exons from Non-infected cells:
# qval_cutoff<-1.30103
# vero_peaks<-fread("/peakCalling/peaks_infectedCellsAll_callsummits_pvalue0.01/gathered_single_peak_files/peaks_noecoli/peaks_nodup_uniq_vero_infectedAll_peaks_all_noecoli.narrowPeak",header = F,sep="\t", col.names = cols) %>% filter(`-log10(q.value)`>qval_cutoff)
# 
# vero_reference_peaks<-vero_peaks %>% mutate(type="retained",condition="vero", non_infected="vero") %>% dplyr::select(type,condition,name,non_infected)

#_________________________________________Metagene plots____________________________________________#
################
#library(metagene2)
inputdir <- "/mappings/covid19_infected_cells/nodup_uniq_alns/all/"
bamfiles<-list.files(inputdir,pattern = ".chroms.nodup.uniq.sorted.bam$")

bamfiles_m6a<- grep("(m6a)",bamfiles,perl = T,ignore.case = T,value = T)
bamfiles_m6a<-bamfiles_m6a[-4]
names(bamfiles_m6a) <-c("WU","SA","UK","Vero") 


regions<- deuGenes_with_peaks_annotation %>% dplyr::select(chr_gtf,start_gtf,end_gtf,groupID,score_gtf,strand_gtf) %>% as.data.frame()
regions.gr<- makeGRangesFromDataFrame(regions,start.field = "start_gtf",end.field = "end_gtf",strand.field = "strand_gtf", seqnames.field = "chr_gtf")
mg <- metagene2$new(regions = regions.gr,bam_files = bamfiles_m6a,  strand_specific = T, invert_strand = T)
mg$produce_metagene()


regions_deu_exons_intersect_m6a<-dex_filt_all_exonic_m6a_intersections %>% dplyr::select(chr_annot,start_annot,end_annot,transcript_id,score_annot,strand_annot)
regions_deu_exons_intersect_m6a.gr<-makeGRangesFromDataFrame(regions_deu_exons_intersect_m6a,start.field = "start_annot",end.field = "end_annot",strand.field = "strand_annot", seqnames.field = "chr_annot")
mg_deu_exons_intersect_m6a<-metagene2$new(regions = regions_deu_exons_intersect_m6a.gr,bam_files = bamfiles_m6a, strand_specific = T, invert_strand = T)
mg_deu_exons_intersect_m6a$produce_metagene()

deuGenes_bed <-flatGTF %>% 
  filter(IDtype_gtf=="aggregate_gene" & features_gtf %in% dex_filt$groupID) %>%
  dplyr::select(chr_gtf,start_gtf,end_gtf,features_gtf,score_gtf,strand_gtf) %>%
  as.data.frame()

deuGenes_regions<-makeGRangesFromDataFrame(deuGenes_bed,start.field = "start_gtf",end.field = "end_gtf",strand.field = "strand_gtf", seqnames.field = "chr_gtf")
mg_deuGenes <-metagene2$new(regions = deuGenes_regions ,bam_files = bamfiles_m6a, strand_specific = T, invert_strand = T)
mg_deuGenes$produce_metagene()

deuSpliceSites_coords<-fread("/plots/manuscript_figures/tables/bedtools_intersect_DEUgenes__splice_sites_coordinates.bed", header = F,sep="\t")
colnames(deuSpliceSites_coords)<-c("chr","start","end","strand")
deuSpliceSites_coords <- deuSpliceSites_coords %>%  mutate(start=start-250,end=end+250)

deuSpliceSites_regions<-makeGRangesFromDataFrame(deuSpliceSites_coords, start.field = "start",end.field = "end",strand.field = "strand")

mg_deuSpliceSites <-metagene2$new(regions = deuSpliceSites_regions ,bam_files = bamfiles_m6a, strand_specific = T, invert_strand = T)

mg_deuSpliceSites$produce_metagene()


############### 
#vennDiagrams: Comparisons at exon level across strains:

deuGeneTxExon_vero <- dex_filt_all_exonic_m6a_intersections %>% filter(condition_peaks=="Vero") %>% dplyr::select(groupID,transcript_id,featureID) %>% mutate(label=paste(groupID,transcript_id,featureID,sep = "_")) %>% dplyr::select(label) %>% distinct() %>% unlist()
names(deuGeneTxExon_vero) <- deuGeneTxExon_vero

deuGeneTxExon_wu <- dex_filt_all_exonic_m6a_intersections %>% filter(condition_peaks=="WU") %>% dplyr::select(groupID,transcript_id,featureID) %>% mutate(label=paste(groupID,transcript_id,featureID,sep = "_")) %>% dplyr::select(label) %>% distinct() %>% unlist()
names(deuGeneTxExon_wu) <- deuGeneTxExon_wu

deuGeneTxExon_uk <- dex_filt_all_exonic_m6a_intersections %>% filter(condition_peaks=="UK") %>% dplyr::select(groupID,transcript_id,featureID) %>% mutate(label=paste(groupID,transcript_id,featureID,sep = "_")) %>% dplyr::select(label) %>% distinct() %>% unlist()
names(deuGeneTxExon_uk) <- deuGeneTxExon_uk

deuGeneTxExon_sa <- dex_filt_all_exonic_m6a_intersections %>% filter(condition_peaks=="SA") %>% dplyr::select(groupID,transcript_id,featureID) %>% mutate(label=paste(groupID,transcript_id,featureID,sep = "_")) %>% dplyr::select(label) %>% distinct() %>% unlist()
names(deuGeneTxExon_sa) <- deuGeneTxExon_sa


ggvenn(list("DEU exons & m6a \n Vero"=deuGeneTxExon_vero,"DEU exons & m6a \n WU"=deuGeneTxExon_wu,"DEU exons & m6a \n UK"=deuGeneTxExon_uk,"DEU exons & m6a \n SA"=deuGeneTxExon_sa))

#heatmap DEU exons with m6a
#dplyr::select(Vero,WU,UK,SA) %>%
dex_filt_all_exonic_m6a_intersect.m <- dex_filt_all_exonic_m6a_intersections %>% 
  dplyr::select(log2fold_WU_Vero,log2fold_UK_Vero,log2fold_SA_Vero) %>% 
  as.matrix()

rownames(dex_filt_all_exonic_m6a_intersect.m)<-paste0(dex_filt_all_exonic_m6a_intersections$gene,"_",dex_filt_all_exonic_m6a_intersections$transcript_id,dex_filt_all_exonic_m6a_intersections$featureID)

pheatmap(dex_filt_all_exonic_m6a_intersect.m, inferno(11), fontsize_row = 2,cluster_cols = F)



# SuppFig4D: Pie Chart DEU genes and m6a levels --------------------------------------
#Number of DEU genes which doesn't have a peak:

vero_m6aDEUGenes<-countsPeaksPerGeneDEU %>% ungroup() %>% filter(condition_peak=="vero") %>% dplyr::select(groupID) %>% distinct() %>% unlist()
deuGenes_nopeak<-setdiff(unique(dex_filt$groupID),unique(vero_m6aDEUGenes))

length(setdiff(unique(dex_filt$groupID),unique(countsPeaksPerGeneDEU$groupID)))



non_m6aDEUGenes<-data.frame(matrix(ncol = 3,nrow = length(deuGenes_nopeak)))
colnames(non_m6aDEUGenes)<-c("condition_peak","groupID","n")
non_m6aDEUGenes$condition_peak <- rep("vero",length(deuGenes_nopeak))
non_m6aDEUGenes$groupID=deuGenes_nopeak
non_m6aDEUGenes$n=rep(0,length(deuGenes_nopeak))

#Bind data frame of non-m6a DEU genes with DEU genes having m6A for categorizing and counting:
mergedCountsPeaksPerGeneDEU <- rbind(as.data.frame(countsPeaksPerGeneDEU),non_m6aDEUGenes)

mergedCountsPeaksPerGeneDEU_summary <- mergedCountsPeaksPerGeneDEU %>% 
  mutate(m6a_level=if_else(n<2 & n>0,"low",
                           if_else(n==2,"medium",
                                   if_else(n>2,"high",
                                           if_else(n==0,"non-m6A","NA"))))) %>% ungroup() %>%
  filter(condition_peak=="vero") %>%  dplyr::select(condition_peak,m6a_level) %>%  dplyr::count(condition_peak,m6a_level)

deu_m6a_proportions <- mergedCountsPeaksPerGeneDEU_summary %>% 
  arrange(desc(m6a_level)) %>%
  mutate(prop = n / sum(mergedCountsPeaksPerGeneDEU_summary$n)) %>%
  mutate(ypos = cumsum(prop)- 0.5*prop )

deu_m6a_proportions <- deu_m6a_proportions %>% mutate(label=paste0(m6a_level,"\n","(",n,")\n",round(prop*100,digits = 2), "%"))

deu_m6a_proportions$m6a_level <- factor(deu_m6a_proportions$m6a_level, levels = c("non-m6A","low","medium","high"))

pdf("/plots/compilation/suppFig4_m6a_levels_in_DEUgenes_pieChart.pdf", width = 8, height = 8)


deu_m6a_proportions %>% 
  ggplot(aes(x="",y=prop*100,fill=m6a_level)) +
  geom_bar(stat="identity", width=1,color="white") +
  coord_polar("y",start = 0) +
  theme_void() +
  theme(legend.position = "none") +
  #geom_label_repel(data=deu_m6a_proportions, aes(label=label, nudge_x=1.5)) +
  geom_text(aes(label=label),color="black",size=4.5, position = position_stack(vjust = 0.5)) +
  scale_fill_manual(values = c("lightgray","#bdd7e7","#6baed6", "#2171b5")) + labs(title = "m6A levels in DEU genes")


dev.off()
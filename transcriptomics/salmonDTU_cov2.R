library(data.table)
library(tidyverse)
library(tximport)
library(rnaseqDTU)
library(GenomicFeatures)
library(DRIMSeq)
library(stageR)
library(DEXSeq)
library(cowplot)


#Prepare annotation (outside R) for DEXSeq differential exon usage analysis:
#python /home/akram/R/x86_64-pc-linux-gnu-library/3.6/DEXSeq/python_scripts/dexseq_prepare_annotation.py chlSab2_wuhCor1_ecoliK12_ncbiGenes_concat.gtf chlSab2_wuhCor1_ecoliK12_ncbiGenes_concat_collapsed_exon_bins.gf

#python /sw/apps/R_packages/4.0.4/rackham/DEXSeq/python_scripts/dexseq_prepare_annotation.py chlSab2_wuhCor1_ecoliK12_ncbiGenes_chroms_concat.gtf chlSab2_wuhCor1_ecoliK12_ncbiGenes_chroms_concat.gff

#/crex/proj/nb_storage/private/BEA21P074_Roshan_Vaid/reference_genomes/Chlorocebus_sabaeus.ChlSab1.1.104.chr_ensembl.gff

dir<-"plots/manuscript_figures/salmon"
#dir="/crex/proj/nb_storage/private/BEA21P074_Roshan_Vaid/salmon_quant"
#samps<-fread("/crex/proj/nb_storage/private/BEA21P074_Roshan_Vaid/salmon_quant/samples.tsv",header = T,sep="\t")
samps<-fread("plots/manuscript_figures/salmon/samples.tsv",header = T,sep="\t")

samps<-as.data.frame(samps)

files<-file.path(dir,samps$sample_id,"quant.sf")
names(files)<-c(paste0(c("Vero","WU","UK","SA"),"_1"),paste0(c("Vero","WU","UK","SA"),"_2"))
files
txi<-tximport(files,type = "salmon",txOut = TRUE,countsFromAbundance = "scaledTPM")

samps$sample_id<-c(paste0(c("Vero","WU","UK","SA"),"_1"),paste0(c("Vero","WU","UK","SA"),"_2"))
samps$condition<-factor(samps$condition)
samps$condition<-relevel(samps$condition,"vero")
#Transcript-to-gene mapping

#gtf<-"reference_genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_concat.gtf"
gtf<-"reference_genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_chroms_concat.gtf"
#gtf="/crex/proj/nb_project/private/genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_concat.gtf"
txdb.filename="reference_genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_concat_txdb.sqlite"
txdb<-makeTxDbFromGFF(gtf,format = "gtf")
saveDb(txdb,txdb.filename)
txdb<-loadDB(txdb.filename)

txdf<-AnnotationDbi::select(txdb,keys(txdb,"GENEID"),"TXNAME","GENEID")
tab<-table(txdf$GENEID)
txdf$ntx<-tab[match(txdf$GENEID,names(tab))]


cov2genes<-c("S","ORF8","ORF7a","ORF7b","ORF6","ORF3a","ORF1ab","ORF10","N","M","E")
cov2_ecoli_genes<-c(cov2genes,grep("gene-",rownames(cts),value=T))

txdf<-txdf %>% filter(!(GENEID %in% cov2_ecoli_genes))

#Filtering by relative abundance:
#txi$abundance %>% as.data.frame() %>%mutate(gene=rownames(.)) %>%  pivot_longer(names_to = "condition",values_to = "abundance",-gene) %>% group_by(gene) %>% mutate(rel_abundance=abundance/sum(abundance)) %>% filter(rel_abundance>0.05) %>% dplyr::select(gene) %>% distinct() %>% dim()

cts<-txi$counts
cts<-cts[which(is.na(match(rownames(cts),cov2_ecoli_genes))),]
cts<-cts[rowSums(cts)>0,]
head(cts)
all(rownames(cts) %in% txdf$TXNAME)

counts<-inner_join(txdf %>% dplyr::select(GENEID,TXNAME),as.data.frame(cts) %>% mutate(TXNAME=rownames(.)),by=c("TXNAME"="TXNAME"))
colnames(counts)<-c("gene_id","feature_id",samps$sample_id)

#Make counts object
d <- dmDSdata(counts=counts, samples=samps)
d

methods(class=class(d))
#Filtering of lowly expressed transcripts 
n<-8
n.small<-2

d<-dmFilter(d,min_samps_feature_expr=n.small,
            min_feature_expr=0.1,
            min_samps_feature_prop=n.small,
            min_feature_prop=0,
            min_samps_gene_expr=n,
            min_gene_expr=1)

d

#Count how many genes have N isoforms:
table(table(counts(d)$gene_id))

#Create a design matrix using the desing formula specified in the samples file:
design_full <- model.matrix(~condition, data=DRIMSeq::samples(d))

colnames(design_full)
#Estimate model parameters for DTU.
#Estimating precision (related with dispersion parameter in Drichlet Multinomial Model), higher dispersion in counts -> less precision.
d <- dmPrecision(d, design=design_full)
d <- dmFit(d, design=design_full)
d <- dmTest(d, coef=c("conditionwu","conditionuk","conditionsa"))
#d <- dmTest(d, coef="condition")

#Statistical analysis of differential transcript usage
res <- DRIMSeq::results(d)
head(res)

res.txp <- DRIMSeq::results(d, level="feature")
head(res.txp)

no.na <- function(x) ifelse(is.na(x), 1, x)
res$pvalue <- no.na(res$pvalue)
res.txp$pvalue <- no.na(res.txp$pvalue)

idx <- which(res$adj_pvalue < 0.05)[1]
res[idx,]

plotProportions(d, res[idx,"gene_id"], "condition")

#stageR screening step, validating evidence of DTU at gene and transcript levels:
pScreen <- res$pvalue
strp <- function(x) substr(x,1,15)
names(pScreen) <- strp(res$gene_id)

pConfirmation <- matrix(res.txp$pvalue, ncol=1)
rownames(pConfirmation) <- strp(res.txp$feature_id)

tx2gene <- res.txp[,c("feature_id", "gene_id")]
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=FALSE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  drim.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                  onlySignificantGenes=TRUE)
})
head(drim.padj)
#The transcripts with values in the column, transcript, less than 0.05 pass the confirmation stage on a target 5% overall false
#discovery rate, or OFDR



res.txp.filt <- DRIMSeq::results(d, level="feature")
smallProportionSD <- function(d, filter=0.1) {
  cts <- as.matrix(subset(counts(d), select=-c(gene_id, feature_id)))
  gene.cts <- rowsum(cts, counts(d)$gene_id)
  total.cts <- gene.cts[match(counts(d)$gene_id, rownames(gene.cts)),]
  props <- cts/total.cts
  propSD <- sqrt(rowVars(props))
  propSD < filter
}
filt <- smallProportionSD(d)
res.txp.filt$pvalue[filt] <- 1 
res.txp.filt$adj_pvalue[filt] <- 1

save(res.txp.filt,file = "/dexseq_counts/res.txp.filt_DTU_DRIMseq.rda")

#DExseq

sample.data <- DRIMSeq::samples(d)
sample.data$condition<-factor(sample.data$condition)
count.data <- round(as.matrix(counts(d)[,-c(1:2)]))
dxd <- DEXSeqDataSet(countData=count.data,
                     sampleData=sample.data,
                     design=~sample + exon + condition:exon,
                     featureID=counts(d)$feature_id,
                     groupID=counts(d)$gene_id)

dxd <- estimateSizeFactors(dxd)
dxd <- estimateDispersions(dxd, quiet=TRUE)
dxd <- testForDEU(dxd, reducedModel=~sample + exon)

dxr <- DEXSeqResults(dxd, independentFiltering=FALSE)
qval <- perGeneQValue(dxr)
dxr.g <- data.frame(gene=names(qval),qval)

columns <- c("featureID","groupID","pvalue")
dxr <- as.data.frame(dxr[,columns])
head(dxr)


strp <- function(x) substr(x,1,15)
pConfirmation <- matrix(dxr$pvalue,ncol=1)
dimnames(pConfirmation) <- list(strp(dxr$featureID),"transcript")
pScreen <- qval
names(pScreen) <- strp(names(pScreen))
tx2gene <- as.data.frame(dxr[,c("featureID", "groupID")])
for (i in 1:2) tx2gene[,i] <- strp(tx2gene[,i])

stageRObj <- stageRTx(pScreen=pScreen, pConfirmation=pConfirmation,
                      pScreenAdjusted=TRUE, tx2gene=tx2gene)
stageRObj <- stageWiseAdjustment(stageRObj, method="dtu", alpha=0.05)
suppressWarnings({
  dex.padj <- getAdjustedPValues(stageRObj, order=FALSE,
                                 onlySignificantGenes=TRUE)
})

txi.g <- tximport(files, type="salmon", tx2gene=txdf[,2:1])
samps$condition<-factor(samps$condition)
samps$condition<-relevel(samps$condition,"vero")
dds <- DESeqDataSetFromTximport(txi.g, samps, ~condition)
dds <- DESeq(dds)
dres <- DESeq2::results(dds)

all(dxr.g$gene %in% rownames(dres))

dres <- dres[dxr.g$gene,]
# we can only color because we simulated...

bigpar()

# here cap the smallest DESeq2 adj p-value
cap.padj <- pmin(-log10(dres$padj), 100)
# this vector only used for plotting
jitter.padj <- -log10(dxr.g$qval + 1e-20)
jp.idx <- jitter.padj == 20
jitter.padj[jp.idx] <- rnorm(sum(jp.idx),20,.25)
plot(cap.padj, jitter.padj, col="blue",
     xlab="Gene expression",
     ylab="Transcript usage")
# legend("topright",
#        c("DGE","DTE","DTU","null"),
#        col=c(1:3,8), pch=20, bty="n")


dxr.sig.genes<- dxr.g %>% filter(qval < 0.01)
dres.drx.sig.genes<- dres[dxr.sig.genes$gene,]
dres.drx.sig.genes

res.txp.filt %>% as.data.frame() %>% filter(adj_pvalue < 0.05) %>% filter(gene_id %in% common_diffexp_genes) -> res.txp.in.commonDEGS
res.txp.in.commonDEGS.genenames<-unique(res.txp.in.commonDEGS$gene_id)

pdf(/plots/manuscript_figures/differentialTrancriptUsage_genes_intersectWith_commonDEGs.pdf",onefile = T)
plotProportions(d,res.txp.in.commonDEGS.genenames[1] , "condition") -> p1
plotProportions(d,res.txp.in.commonDEGS.genenames[2] , "condition") -> p2
plotProportions(d,res.txp.in.commonDEGS.genenames[3] , "condition") -> p3
plotProportions(d,res.txp.in.commonDEGS.genenames[4] , "condition") -> p4
plotProportions(d,res.txp.in.commonDEGS.genenames[5] , "condition") -> p5
plotProportions(d,res.txp.in.commonDEGS.genenames[6] , "condition") -> p6
plotProportions(d,res.txp.in.commonDEGS.genenames[7] , "condition") -> p7
plotProportions(d,res.txp.in.commonDEGS.genenames[8] , "condition") -> p8
p1
p2
p3
p4
p5
p6
p7
p8


dev.off()

#ggsave(differentialTrancriptUsage_genes_intersectWith_commonDEGs.pdf", device = "pdf", width = 8, height = 6,units = "in")


#Save drim padj datasset with filtered genes that passed DTU testing with DRIMseq
openxlsx::write.xlsx(list("DRIMseq_padj"=res.txp.filt %>% filter(adj_pvalue < 0.05)),/plots/manuscript_figures/tables/DTU_drimseq_results.xls")

res.txp.significant.genes <- res.txp.filt %>% filter(adj_pvalue < 0.05) %>% dplyr::select(gene_id) %>% distinct() %>% unlist()

names(res.txp.significant.genes)<-res.txp.significant.genes

res.txp.significant.genes_entrez <-mapIds(org.Hs.eg.db,keys = res.txp.significant.genes, column = "ENTREZID",keytype = "SYMBOL",multiVals = "first")
names(res.txp.significant.genes_entrez)<-res.txp.significant.genes_entrez
# Enrichment analysis
#####

goBP_res.txp.significant.genes<-enrichGO(gene=res.txp.significant.genes,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="BP",pvalueCutoff = 0.05,readable = F)%>% gofilter(.,level = c(7:10))
goMF_res.txp.significant.genes<-enrichGO(gene=res.txp.significant.genes,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="MF",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
goCC_res.txp.significant.genes<-enrichGO(gene=res.txp.significant.genes,OrgDb = org.Hs.eg.db, keyType = "SYMBOL", ont="CC",pvalueCutoff = 0.05,readable = F) %>% gofilter(.,level = c(7:10))
kegg_res.txp.significant.genes<-enrichKEGG(res.txp.significant.genes_entrez,pvalueCutoff = 0.05)
reactome_res.txp.significant.genes<-enrichPathway(res.txp.significant.genes_entrez,pvalueCutoff=0.05)
enrichrWP_res.txp.significant.genes<-enrichr(res.txp.significant.genes,"WikiPathways_2019_Human")
enrichrCov2sets_res.txp.significant.genes<-enrichr(res.txp.significant.genes,"COVID-19_Related_Gene_Sets")

pdf(/plots/manuscript_figures/enrichmentAnalysis_genesWith_DiffTranscriptUsage.pdf", width = 8, height = 16)
mydotplot(goBP_res.txp.significant.genes,showCategory=30,x="Count", title="GO:Biological Process \n Genes with differential transcript usage")
mydotplot(goMF_res.txp.significant.genes,showCategory=30,x="Count",title="GO:Molecular function \n Genes with differential transcript usage")
mydotplot(goCC_res.txp.significant.genes,showCategory=30,x="Count", title="GO:Cellular component \n Genes with differential transcript usage")
mydotplot(kegg_res.txp.significant.genes,showCategory=30,x="Count", title="KEGG pathways \n Genes with differential transcript usage")

plotEnrich(enrichrCov2sets_res.txp.significant.genes$`COVID-19_Related_Gene_Sets` %>% filter(P.value<0.05)) + labs(title="Enrichr: COVID-19 Related Gene Sets \ Genes with differential transcript usage ")
plotEnrich(enrichrWP_res.txp.significant.genes$WikiPathways_2019_Human %>% filter(P.value<0.05))  + labs(title="Enrichr: WikiPathways \ Genes with differential transcript usage ")
dev.off()

goTables_dtu<-list("Biological_Process"=goBP_res.txp.significant.genes@result %>% filter(pvalue < 0.05),
                   "Molecular function"=goMF_res.txp.significant.genes@result %>% filter(pvalue < 0.05),
                   "Cellular component"=goCC_res.txp.significant.genes@result %>% filter(pvalue < 0.05),
                   "KEGG"=kegg_res.txp.significant.genes@result %>% filter(pvalue < 0.05),
                   "Enrichr_Cov2sets"=enrichrCov2sets_res.txp.significant.genes$`COVID-19_Related_Gene_Sets` %>% filter(P.value < 0.05),
                   "Enrichr_WikiPathways"=enrichrWP_res.txp.significant.genes$WikiPathways_2019_Human %>% filter(P.value < 0.05))
openxlsx::write.xlsx(goTables_dtu,/plots/manuscript_figures/tables/enrcihmentAnalysis_differentialTranscriptUsage_all.xls")


#Exploring proportions dataset
##### 

proportions(d) %>% filter(gene_id %in% res.txp.in.commonDEGS.genenames) %>% gather(key="sample",value="proportion",-gene_id,-feature_id) %>% mutate(gene=paste(gene_id,feature_id,sep="\n"),condition=gsub("(_1|_2)","",sample,perl = T)) %>% ggplot(aes(gene,proportion,fill=condition)) + geom_bar(stat = "identity",position = position_stack()) + coord_flip()

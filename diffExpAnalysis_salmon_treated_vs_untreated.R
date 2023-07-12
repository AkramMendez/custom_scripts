library(DESeq2)
library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(data.table)
library(ggtext)
library(factoextra)
library(pheatmap)
library(cowplot)
library(EnhancedVolcano)
library(sva)
library(tximport)
library(tximeta)
library(tidyverse)
library(GenomicFeatures)
library(ensembldb)
library(EnsDb.Hsapiens.v86)
library(rtracklayer)
library(org.Hs.eg.db)

# Import salmon quant files -----------------------------------------------
salmonDir <- "../salmon_quant"
salmon_samples <- list.files(salmonDir,recursive = T,pattern = "quant.sf",full.names = T)

# Create tx2gene table ----------------------------------------------------
# Set up filenames. Only change this line to point to the GTF.
# Use the same annotation file as the one used for Salmon quantification:
gtf <-"gencode.v38.annotation.gff3"
#"gencode.v36.chr_patch_hapl_scaff.basic.annotation.gtf"

# These filenames automatically generated from the above.
txdb.filename <- gsub("gff3", "sqlite", gtf)
tx2gene.filename <- gsub("gff3", "tx2gene.csv", gtf)


# Check to see if the txdb database file exists.
if (!file.exists(txdb.filename)) {
  ## If not, make it. Do it once. This takes a while.
  message(paste(txdb.filename, "doesn't exist. Creating and saving to disk..."))
  txdb <- makeTxDbFromGFF(gtf)
  #txdb <- makeTxDbFromGRanges(gtf)
  saveDb(txdb, txdb.filename)
  message("Done.")
} else {
  ## If it already exists, load it from file, quickly!
  message(paste(txdb.filename, "found. Loading..."))
  txdb <- loadDb(txdb.filename)
  message("Done.")
}

keytypes(txdb)



# Create the tx2gene
tx2gene <- mapIds(txdb,
                  keys=keys(txdb,"GENEID"),
                  column="TXNAME",
                  keytype="GENEID",
                  multiVals="list") %>%
  enframe(name="ensgene", value="enstxp") %>%
  unnest(cols = c(enstxp)) %>%
  dplyr::select(enstxp, ensgene) %>%
  distinct()

#Translate Ensemble IDs to gene and transcript names names:


convertEnsemblToGene <- function(ensemblIds){
  ensemblIdstrim <- unique(gsub("\\.\\d+","", ensemblIds))
  #(run once)
  library(biomaRt)
  mart <- useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
  genes <-getBM(filters = "ensembl_gene_id", attributes = c("ensembl_gene_id","hgnc_symbol"),ensemblIdstrim,mart = mart)
  #genes <- mapIds(org.Hs.eg.db, keys = ensemblIdstrim,keytype = "ENSEMBL", column = "SYMBOL", multiVals = "first")
  #df <- data.frame(ensembl_gene=ensemblIdstrim, gene=genes)
  
  return(genes)
}

# Generate table of EnsembleIDs and gene symbols
ensemblToGene <- convertEnsemblToGene(tx2gene$ensgene)

#_________________________________________________________________________________________________________________________________________________________________________________#
#_____________________________________________________________________________________#
# Differential expression Treat vs Control:
#_____________________________________________________________________________________#

# Import Salmon quantifications -------------------------------------------

coldata_cellType1 <- fread("../salmonSamplesInfo_coldata.csv", sep=",", header = T)


# Import salmon counts
salmon_table_cellType1 <- tximeta(coldata_cellType1,type = "salmon", skipSeqinfo = T)
salmon_table_cellType1 <- addIds(salmon_table_cellType1,"SYMBOL",gene = T)
geneLevelSalmon_cellType1 <- summarizeToGene(salmon_table_cellType1)


ddsSalmon_cellType1_untreated_kdvsctrl <- DESeqDataSet(geneLevelSalmon_cellType1, design = ~ type)

# Assign levels to reference conditions:
ddsSalmon_cellType1_untreated_kdvsctrl$type <- factor(ddsSalmon_cellType1_untreated_kdvsctrl$type, c("controlsh","kd"))

ddsSalmon_cellType1_untreated_kdvsctrl <- DESeq(ddsSalmon_cellType1_untreated_kdvsctrl,modelMatrixType = "standard", fitType = "local")

resultsNames(ddsSalmon_cellType1_untreated_kdvsctrl)

results(ddsSalmon_cellType1_untreated_kdvsctrl) %>% as.data.frame() %>% arrange(log2FoldChange)


normCounts_kd_vs_ctrl_untreated_cellType1 <- counts(ddsSalmon_cellType1_untreated_kdvsctrl, normalized=T)

normCounts_kd_vs_ctrl_untreated_cellType1 <- left_join(
  normCounts_kd_vs_ctrl_untreated_cellType1 %>% as.data.frame() %>% mutate(ensembl=gsub("\\.\\d+","",rownames(.),perl = T)),
  ensemblToGene, by=c("ensembl"="ensembl_gene_id")) %>% mutate(hgnc_symbol=if_else(is.na(hgnc_symbol) | hgnc_symbol =="", ensembl,hgnc_symbol))

left_join(
  normCounts_kd_vs_ctrl_untreated_cellType1 %>% as.data.frame() %>% mutate(ensembl=gsub("\\.\\d+","",rownames(.),perl = T)),
          ensemblToGene, by=c("ensembl"="ensembl_gene_id")) %>% mutate(hgnc_symbol=if_else(is.na(hgnc_symbol) | hgnc_symbol =="", ensembl,hgnc_symbol)) %>% mutate(description="NA") %>% dplyr::select(hgnc_symbol,description,starts_with("kd"),starts_with("control")) %>% dplyr::filter(!grepl("ENSG\\d+",hgnc_symbol, perl = T)) %>% dplyr::filter(!duplicated(hgnc_symbol)) %>% fwrite(., "~/MondalLab/NCC_RNAseq/diffExpAnalysis/tables/be2_tncc_batch/GSEA/normCounts_kd_vs_ctrl_untreated_cellType1.txt", sep = "\t",col.names = T)


lfcShrinkSalmon_cellType1_untreated_kdvsctrl <- lfcShrink(ddsSalmon_cellType1_untreated_kdvsctrl,coef = "type_kd_vs_controlsh",type = 'apeglm') %>% as.data.frame() %>%  mutate(ensembl_gene=gsub("\\.\\d+","",rownames(.),perl = T)) %>% mutate(contrast="kd_untreated_vs_controlsh_untreated")

# Add gene symbols:

lfcShrinkSalmon_cellType1_untreated_kdvsctrl <- left_join(lfcShrinkSalmon_cellType1_untreated_kdvsctrl, ensemblToGene, by=c("ensembl_gene"="ensembl_gene_id"))

lfcShrinkSalmon_cellType1_untreated_kdvsctrl <- lfcShrinkSalmon_cellType1_untreated_kdvsctrl %>% mutate(gene=if_else(hgnc_symbol=="" | is.na(hgnc_symbol),ensembl_gene,hgnc_symbol)) %>% dplyr::select(-hgnc_symbol)


# Check counts
plotCounts(ddsSalmon_cellType1_untreated_kdvsctrl,gene = "ENSG00000165819.12", returnData=F,intgroup = "type", normalized = T)  
plotCounts(ddsSalmon_cellType1_untreated_kdvsctrl,gene = "ENSG00000134323.12", returnData=F,intgroup = "type", normalized = T)


#####

#_________________________________________________________________________________________________________________________________________________________________________________#
# Differential expression KD vs Ctrl in treated cells:
#_____________________________________________________________________________________#


coldata_cellType1_Treat <- fread("../salmonSamplesInfo_coldata.csv", sep=",", header = T)


# Import salmon counts for CellType1 samples:
salmon_table_cellType1_Treat <- tximeta(coldata_cellType1_Treat,type = "salmon", skipSeqinfo = T)
salmon_table_cellType1_Treat <- addIds(salmon_table_cellType1_Treat,"SYMBOL",gene = T)
geneLevelSalmon_cellType1_Treat <- summarizeToGene(salmon_table_cellType1_Treat)


ddsSalmon_cellType1_Treat_kdvsctrl <- DESeqDataSet(geneLevelSalmon_cellType1_Treat, design = ~ type)

# Assign levels to reference conditions:
ddsSalmon_cellType1_Treat_kdvsctrl$type <- factor(ddsSalmon_cellType1_Treat_kdvsctrl$type, c("controlsh","kd"))

ddsSalmon_cellType1_Treat_kdvsctrl <- DESeq(ddsSalmon_cellType1_Treat_kdvsctrl,modelMatrixType = "standard", fitType = "local")

resultsNames(ddsSalmon_cellType1_Treat_kdvsctrl)

results(ddsSalmon_cellType1_Treat_kdvsctrl) %>% as.data.frame() %>% arrange(log2FoldChange)


lfcShrinkSalmon_cellType1_Treat_kdvsctrl <- lfcShrink(ddsSalmon_cellType1_Treat_kdvsctrl,coef = "type_kd_vs_controlsh",type = 'apeglm') %>% as.data.frame() %>%  mutate(ensembl_gene=gsub("\\.\\d+","",rownames(.),perl = T)) %>% mutate(contrast="kd_Treattreated_vs_controlsh_Treattreated")

# Add gene symbols:

lfcShrinkSalmon_cellType1_Treat_kdvsctrl <- left_join(lfcShrinkSalmon_cellType1_Treat_kdvsctrl, ensemblToGene, by=c("ensembl_gene"="ensembl_gene_id"))

lfcShrinkSalmon_cellType1_Treat_kdvsctrl <- lfcShrinkSalmon_cellType1_Treat_kdvsctrl %>% mutate(gene=if_else(hgnc_symbol=="" | is.na(hgnc_symbol),ensembl_gene,hgnc_symbol)) %>% dplyr::select(-hgnc_symbol)


lfcShrinkSalmon_cellType1_Treat_kdvsctrl %>% fwrite(.,"lfcShrinkSalmon_cellType1_Treat_kdvsctrl.tsv", sep="\t")

# Check kd counts
plotCounts(ddsSalmon_cellType1_Treat_kdvsctrl,gene = "ENSG00000165819.12", returnData=F,intgroup = "type", normalized = T) 

plotCounts(ddsSalmon_cellType1_Treat_kdvsctrl,gene = "ENSG00000134323.12", returnData=F,intgroup = "type", normalized = T) 

#_________________________________________________________________________________________________________________________________________________________________________________#
# Differential expression 
# Control and KD in
# Treated vs Untreated CellType1 cells
#_____________________________________________________________________________________#

###### ------------------------- Treated vs Untreated Control CellType1 ---------------------------------######

coldata_cellType1_TreattreatedVSuntreated <- fread("salmonSamplesInfo_coldata_cellType1_Treated.csv", sep=",", header = T)


# Import salmon counts for CellType1 samples:
salmon_table_cellType1_control_TreattreatedVSuntreated <- tximeta(
  coldata_cellType1_TreattreatedVSuntreated %>% filter(type == "controlsh"),
  type = "salmon",
  skipSeqinfo = T
)
salmon_table_cellType1_control_TreattreatedVSuntreated <- addIds(salmon_table_cellType1_control_TreattreatedVSuntreated,"SYMBOL",gene = T)
geneLevelSalmon_cellType1_control_TreattreatedVSuntreated <- summarizeToGene(salmon_table_cellType1_control_TreattreatedVSuntreated)

# Differential expression analysis 
ddsSalmon_cellType1_control_TreattreatedVSuntreated <- DESeqDataSet(geneLevelSalmon_cellType1_control_TreattreatedVSuntreated, design = ~ condition) #~ type + condition

# Assign levels to reference conditions:
ddsSalmon_cellType1_control_TreattreatedVSuntreated$type <- factor(ddsSalmon_cellType1_control_TreattreatedVSuntreated$type, c("controlsh","kd"))
ddsSalmon_cellType1_control_TreattreatedVSuntreated$condition <- factor(ddsSalmon_cellType1_control_TreattreatedVSuntreated$condition, c("untreated","treated"))

ddsSalmon_cellType1_control_TreattreatedVSuntreated <- DESeq(ddsSalmon_cellType1_control_TreattreatedVSuntreated,modelMatrixType = "standard", fitType = "local")

resultsNames(ddsSalmon_cellType1_control_TreattreatedVSuntreated)

results(ddsSalmon_cellType1_control_TreattreatedVSuntreated) %>% as.data.frame() %>% arrange(log2FoldChange)


lfcShrinkSalmon_cellType1_control_TreattreatedVSuntreated <- lfcShrink(ddsSalmon_cellType1_control_TreattreatedVSuntreated,coef = "condition_Treat_vs_untreated",type = 'apeglm') %>% as.data.frame() %>%  mutate(ensembl_gene=gsub("\\.\\d+","",rownames(.),perl = T)) %>% mutate(contrast="RAtreated_controlsh_VS_untreated_controlsh")

# Add gene symbols:

lfcShrinkSalmon_cellType1_control_TreattreatedVSuntreated <- left_join(lfcShrinkSalmon_cellType1_control_TreattreatedVSuntreated, ensemblToGene, by=c("ensembl_gene"="ensembl_gene_id"))

lfcShrinkSalmon_cellType1_control_TreattreatedVSuntreated <- lfcShrinkSalmon_cellType1_control_TreattreatedVSuntreated %>% mutate(gene=if_else(hgnc_symbol=="" | is.na(hgnc_symbol),ensembl_gene,hgnc_symbol)) %>% dplyr::select(-hgnc_symbol)


# Check counts for single genes:
plotCounts(ddsSalmon_cellType1_control_TreattreatedVSuntreated,gene = "ENSG00000165819.12", returnData=F,intgroup = "condition", normalized = T)
plotCounts(ddsSalmon_cellType1_control_TreattreatedVSuntreated,gene = "ENSG00000134323.12", returnData=F,intgroup = "condition", normalized = T)

plotMA(ddsSalmon_cellType1_control_TreattreatedVSuntreated)


###### ------------------------- Treated vs Untreated KD CellType1 ---------------------------------######


# Import salmon counts for CellType1 samples:
salmon_table_cellType1_kd_TreattreatedVSuntreated <- tximeta(
  coldata_cellType1_TreattreatedVSuntreated %>% filter(type == "kd"),
  type = "salmon",
  skipSeqinfo = T
)

salmon_table_cellType1_kd_TreattreatedVSuntreated <- addIds(salmon_table_cellType1_kd_TreattreatedVSuntreated,"SYMBOL",gene = T)
geneLevelSalmon_cellType1_kd_TreattreatedVSuntreated <- summarizeToGene(salmon_table_cellType1_kd_TreattreatedVSuntreated)

# Differential expression analysis 
ddsSalmon_cellType1_kd_TreattreatedVSuntreated <- DESeqDataSet(geneLevelSalmon_cellType1_kd_TreattreatedVSuntreated, design = ~ condition) #~ type + condition

# Assign levels to reference conditions:
ddsSalmon_cellType1_kd_TreattreatedVSuntreated$type <- factor(ddsSalmon_cellType1_kd_TreattreatedVSuntreated$type, c("controlsh","kd"))
ddsSalmon_cellType1_kd_TreattreatedVSuntreated$condition <- factor(ddsSalmon_cellType1_kd_TreattreatedVSuntreated$condition, c("untreated","treated"))

ddsSalmon_cellType1_kd_TreattreatedVSuntreated <- DESeq(ddsSalmon_cellType1_kd_TreattreatedVSuntreated,modelMatrixType = "standard", fitType = "local")

resultsNames(ddsSalmon_cellType1_kd_TreattreatedVSuntreated)

results(ddsSalmon_cellType1_kd_TreattreatedVSuntreated) %>% as.data.frame() %>% arrange(log2FoldChange)


lfcShrinkSalmon_cellType1_kd_TreattreatedVSuntreated <- lfcShrink(ddsSalmon_cellType1_kd_TreattreatedVSuntreated,coef = "condition_Treat_vs_untreated",type = 'apeglm') %>% as.data.frame() %>%  mutate(ensembl_gene=gsub("\\.\\d+","",rownames(.),perl = T)) %>% mutate(contrast="RAtreated_kd_VS_untreated_kd")

# Add gene symbols:

lfcShrinkSalmon_cellType1_kd_TreattreatedVSuntreated <- left_join(lfcShrinkSalmon_cellType1_kd_TreattreatedVSuntreated, ensemblToGene, by=c("ensembl_gene"="ensembl_gene_id"))

lfcShrinkSalmon_cellType1_kd_TreattreatedVSuntreated <- lfcShrinkSalmon_cellType1_kd_TreattreatedVSuntreated %>% mutate(gene=if_else(hgnc_symbol=="" | is.na(hgnc_symbol),ensembl_gene,hgnc_symbol)) %>% dplyr::select(-hgnc_symbol)


# Check single gene counts 
plotCounts(ddsSalmon_cellType1_kd_TreattreatedVSuntreated,gene = "ENSG00000165819.12", returnData=F,intgroup = "condition", normalized = T)
plotCounts(ddsSalmon_cellType1_kd_TreattreatedVSuntreated,gene = "ENSG00000134323.12", returnData=F,intgroup = "condition", normalized = T)

plotMA(ddsSalmon_cellType1_kd_TreattreatedVSuntreated)

# Save LFC tables:

whichCols <- c("gene","baseMEan","log2FoldChange","lfcSE","pvalue","padj")

list(
  "KDvsCtrl_untreated_CellType1"=signif_cellType1_untreated_kdvsctrl,
  "KDvsCtrl_treated_CellType1"=signif_cellType1_Treat_kdvsctrl
) %>% openxlsx::write.xlsx("diffExpTables_cellType1_KDcsCtrl_treated_untreated.xlsx", overwrite = T)

library(pheatmap)
library(corrplot)
library(sva)

filesPub <- c(
#vero8h_mock_R1="/public_data/salmon_quant/SRR13091744_vero8h_mock_rep1/quant.sf",
#vero8h_mock_R2="/public_data/salmon_quant/SRR13091745_vero8h_mock_rep2/quant.sf",
#vero8h_cov2_R1="/public_data/salmon_quant/SRR13091741_vero8h_cov2Infected_rep1/quant.sf",
#vero8h_cov2_R2="/public_data/salmon_quant/SRR13091742_vero8h_cov2Infected_rep2/quant.sf",
vero24h_mock_R1="/public_data/salmon_quant/SRR12164495_vero24h_mock_rep1/quant.sf",
vero24h_mock_R2="/public_data/salmon_quant/SRR12164495_vero24h_mock_rep2/quant.sf",
vero24h_cov2_R1="/public_data/salmon_quant/SRR12164499_vero24h_cov2Infected_rep2/quant.sf",
vero24h_cov2_R2="/public_data/salmon_quant/SRR12164499_vero24h_cov2Infected_rep2/quant.sf"
)

#vero24h_/tables/InfectedvsNonInfected_veroE6_GSE153940_enrichR.txt

#Vero 8h: PRJNA679776, Zhang Y et al., "SARS-CoV-2 hijacks folate and one-carbon metabolism for viral replication.", Nat Commun, 2021 Mar 15;12(1):1676
#Vero 24h: GSE153940, PRJNA644588, Riva L et al., "Discovery of SARS-CoV-2 antiviral drugs through large-scale compound repurposing.", Nature, 2020 Oct;586(7827):113-119

# Importing counts B.1, B.1.1.7, B.1.1351
dir<-"/plots/manuscript_figures/salmon"
samps<-fread("/plots/manuscript_figures/salmon/samples.tsv",header = T,sep="\t")

samps<-as.data.frame(samps)

files<-file.path(dir,samps$sample_id,"quant.sf")
names(files)<-c(paste0(c("Vero","B.1","B.1.1.7","B.1.351"),"_1"),paste0(c("Vero","B.1","B.1.1.7","B.1.351"),"_2"))
files

txi<-tximport(files,type = "salmon",txOut = TRUE,countsFromAbundance = "scaledTPM")

samps$sample_id<-c(paste0(c("Vero","B.1","B.1.1.7","B.1.351"),"_1"),paste0(c("Vero","B.1","B.1.1.7","B.1.351"),"_2"))
samps$condition<-factor(samps$condition)
samps$condition<-relevel(samps$condition,"vero")

#Transcript-to-gene mapping
#txdb.filename="/reference_genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_concat_txdb.sqlite"

txdb<-loadDb(txdb.filename)

txdf<-AnnotationDbi::select(txdb,keys(txdb,"GENEID"),"TXNAME","GENEID")
tab<-table(txdf$GENEID)
txdf$ntx<-tab[match(txdf$GENEID,names(tab))]


cov2genes<-c("S","ORF8","ORF7a","ORF7b","ORF6","ORF3a","ORF1ab","ORF10","N","M","E", "ORF1a")
cov2_ecoli_genes<-c(cov2genes,grep("gene-",rownames(cts),value=T))

txdf<-txdf %>% filter(!(GENEID %in% cov2_ecoli_genes))

#Filtering by relative abundance:
#txi$abundance %>% as.data.frame() %>%mutate(gene=rownames(.)) %>%  pivot_longer(names_to = "condition",values_to = "abundance",-gene) %>% group_by(gene) %>% mutate(rel_abundance=abundance/sum(abundance)) %>% filter(rel_abundance>0.05) %>% dplyr::select(gene) %>% distinct() %>% dim()

#Import scaledTPM values for correlation matrix or raw counts for differential expression analysis:
#txiPub<-tximport(filesPub,type = "salmon",countsFromAbundance = "scaledTPM", tx2gene = txdf[,c(2,1)], txIn = T)
txiPub <- tximport(filesPub,type = "salmon",countsFromAbundance = "no", tx2gene = txdf[,c(2,1)], txIn = T)

# Merge lab's data with public data (Include only SARS-CoV-2 infected samples:
allData <- inner_join(as.data.frame(txi$abundance) %>% dplyr::select(matches("B.")) %>% rownames_to_column("gene"),
                      as.data.frame(txiPub$abundance) %>% rownames_to_column("gene"), 
                      by=c("gene"="gene")
                      )
batch <- c(rep(1,6),2,2)

adjustBatchAllData <- ComBat_seq(as.matrix(allData[,-1]), batch = batch,group = NULL)

rownames(adjustBatchAllData) <- allData[,1]
colnames(adjustBatchAllData) <- str_wrap(c("Infected-B.1 R1","Infected-B.1.1.7 R1","Infected-B.1.351 R1",
                                           "Infected-B.1 R2","Infected-B.1.1.7 R2","Infected-B.1.351 R2",
                                           "SARS-CoV-2 Riva L et al. 2020 R1", "SARS-CoV-2 Riva L et al. 2020 R2"
),20)

# colnames(adjustBatchAllData) <- str_wrap(c("Infected-B.1 R1","Infected-B.1.1.7 R1","Infected-B.1.351 R1",
#                                   "Infected-B.1 R2","Infected-B.1.1.7 R2","Infected-B.1.351 R2",
#                                   "SARS-CoV-2 Zhang Y et al. 2021 R1", "SARS-CoV-2 Zhang Y et al. 2021 R2",
#                                   "SARS-CoV-2 Riva L et al. 2020 R1", "SARS-CoV-2 Riva L et al. 2020 R2"
#                                   ),20)

#corrAllData <- cor(as.matrix(allData[,-1]))
#corrAllData <- cor(adjustBatchAllData[,c(1:6)])
corrAllData <- cor(adjustBatchAllData)

svg("/figures/correlationSamples_labAndPublicData.svg", width = 12,height = 12)

corrplot(corrAllData,order = 'hclust', 
         addrect = 0, 
         addCoef.col = 'white',
         method = 'circle',
         type="upper",
         number.cex = 0.6,
         tl.srt = 45,
#         col = rev(coolwarm(100)),
         diag = T,
is.corr = T,tl.col = "black")

#colorlegend(xlim=c(10,15), ylim=c(10,15), rev(coolwarm(100)), c(seq(-1,1,.25)), align="l", vertical=TRUE, addlabels=TRUE)

dev.off()

# Sample distances heatmap
sampleDistsPub <- dist(t(adjustBatchAllData))
sampleDistsPub.m <- as.matrix(sampleDistsPub)

pheatmap(sampleDistsPub.m, 
         clustering_distance_rows = sampleDistsPub, 
         clustering_distance_cols = sampleDistsPub,
         col=viridis::viridis(100,option = "G",begin = 0.2,end = 0.8, direction = -1), 
         border_color = NA,
         width = 12, 
         height = 12,
         treeheight_row = 100,
         treeheight_col = 100,
         legend_labels = seq(0,350000,by=100000),
         legend_breaks = seq(0,350000,by=100000)#,
         #filename = "/figures/sampleDistances_labAndPublicData_TPM_batchAdjusted.tiff",
        )

# PCA
fviz_pca_ind(PCA(t(adjustBatchAllData)))



# Differential expression analysis ----------------------------------------

metadata_vero8h <- NULL
metadata_vero8h$sample<-as.factor(c("vero8h_mock_R1","vero8h_mock_R2","vero8h_cov2_R1","vero8h_cov2_R2"))
metadata_vero8h$condition<-as.factor(c("Non_infected","Non_infected","SARSCoV2_infected","SARSCoV2_infected"))
metadata_vero8h$condition<-relevel(metadata_vero8h$condition,"Non_infected")
metadata_vero8h<-as.data.frame(metadata_vero8h)

# metadata_vero24h <- NULL
# metadata_vero24h$sample<-as.factor(c("vero24h_mock_R1","vero24h_mock_R2","vero24h_cov2_R1","vero24h_cov2_R2"))
# metadata_vero24h$condition<-as.factor(c("Non_infected","Non_infected","SARSCoV2_infected","SARSCoV2_infected"))
# metadata_vero24h$condition<-relevel(metadata_vero24h$condition,"Non_infected")
# metadata_vero24h<-as.data.frame(metadata_vero24h)

metadata_vero24h <- NULL
metadata_vero24h$sample<-as.factor(c("vero24h_mock_R1","vero24h_mock_R2","vero24h_cov2_R1","vero24h_cov2_R2"))
metadata_vero24h$condition<-as.factor(c("Non_infected","Non_infected","SARSCoV2_infected","SARSCoV2_infected"))
metadata_vero24h$condition<-relevel(metadata_vero24h$condition,"Non_infected")
metadata_vero24h<-as.data.frame(metadata_vero24h)


#geneIdx <- match(log2FC.all.df$gene ,rownames(txiPub$counts))

# dds_vero8h <- DESeqDataSetFromMatrix(countData = as.data.frame(txiPub$counts) %>% 
#                                        filter(rownames(.) %in% log2FC.all.df$gene) %>% 
#                                        dplyr::select(c(1:4)) %>% as.matrix() %>% round(), 
#                                      colData=metadata_vero8h, design= ~ condition)

dds_vero24h <- DESeqDataSetFromMatrix(countData = as.data.frame(txiPub$counts) %>% 
                                        filter(rownames(.) %in% log2FC.all.df$gene) %>% 
                                        as.matrix() %>% round(), 
                                      colData=metadata_vero8h, design= ~ condition)
#design(dds)<- ~ type

#dds_vero8h<-DESeq(dds_vero8h,modelMatrixType = "standard")

dds_vero24h<-DESeq(dds_vero24h,modelMatrixType = "standard")

#resultsNames(dds_vero8h)

resultsNames(dds_vero24h)


#_________________________________#
#normalizedCounts_vero8h <- counts(dds_vero8h, normalized=TRUE)

#normalizedCounts_vero8h <-normalizedCounts_vero8h %>% as.data.frame() %>%  mutate(gene=rownames(normalized_counts))


diffExp_vero8h <- results(dds_vero8h,contrast = c("condition","SARSCoV2_infected","Non_infected")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="SARSCoV2_infected_vs_Non_infected")

diffExp_vero24h <- results(dds_vero24h,contrast = c("condition","SARSCoV2_infected","Non_infected")) %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="SARSCoV2_infected_vs_Non_infected")

#Shrink LFC estimates:
diffExpLFCshrink_vero8h <- lfcShrink(dds_vero8h,coef = "condition_SARSCoV2_infected_vs_Non_infected",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="condition_SARSCoV2_infected_vs_Non_infected")

diffExpLFCshrink_vero24h <- lfcShrink(dds_vero24h,coef = "condition_SARSCoV2_infected_vs_Non_infected",type = 'apeglm') %>% as.data.frame() %>%  mutate(gene=rownames(.)) %>% mutate(contrast="condition_SARSCoV2_infected_vs_Non_infected")

# significantly changing genes:

signif_vero8h <- diffExpLFCshrink_vero8h %>% filter(padj < 0.1 & !(gene %in% cov2genes)) %>%  arrange(padj) %>% dplyr::select(gene) %>% distinct() %>% unlist()

signif_vero24h <- diffExpLFCshrink_vero24h %>% filter(padj < 0.1 & !(gene %in% cov2genes)) %>%  arrange(padj) %>% dplyr::select(gene) %>% distinct() %>% unlist()

#---------scatter plots log2 fold changes

pairedCor2 <- function(df, labelx,labely, cor.df){
  colnames(df) <- c("gene","x","y")
  df %>% ggplot(aes(x,y)) + 
    geom_point(size=0.01,alpha=0.3) + 
    geom_smooth(method = lm,fill="blue") +
    stat_pvalue_manual(data = cor.df,
                       label = "{cor}{p.adj}",
                       y.position = "y.position", size=3.5) +
    #ylim(c(-2,2)) + xlim(c(-2,2)) +
    coord_cartesian(xlim = c(-5,5),ylim = c(-5,5)) +
    theme(aspect.ratio = 1) +
    labs(x=paste("log2 Fold Change",labelx, sep="\n"),y=paste("log2 Fold Change",labely, sep = "\n")) +
    theme_prism(base_size = 10,border = T) -> p
  
  return(p)
}



lfc_B.1vsVero8h <- inner_join(log2FC.all.df %>% dplyr::select(gene,B.1),
           diffExpLFCshrink_vero8h %>% 
             dplyr::select(gene,log2FoldChange) %>% 
             dplyr::rename("SARS-CoV-2 infected Vero 8h"=log2FoldChange),
           by="gene") 

lfc_B.1.1.7vsVero8h <- inner_join(log2FC.all.df %>% dplyr::select(gene,B.1.1.7),
                              diffExpLFCshrink_vero8h %>% 
                                dplyr::select(gene,log2FoldChange) %>% 
                                dplyr::rename("SARS-CoV-2 infected Vero 8h"=log2FoldChange),
                              by="gene")

lfc_B.1.351vsVero8h <- inner_join(log2FC.all.df %>% dplyr::select(gene,B.1.351),
                                  diffExpLFCshrink_vero8h %>% 
                                    dplyr::select(gene,log2FoldChange) %>% 
                                    dplyr::rename("SARS-CoV-2 infected Vero 8h"=log2FoldChange),
                                  by="gene")

lfc_B.1vsVero24h <- inner_join(log2FC.all.df %>% dplyr::select(gene,B.1),
                              diffExpLFCshrink_vero24h %>% 
                                dplyr::select(gene,log2FoldChange) %>% 
                                dplyr::rename("SARS-CoV-2 infected Vero 24h"=log2FoldChange),
                              by="gene")

lfc_B.1.1.7vsVero24h <- inner_join(log2FC.all.df %>% dplyr::select(gene,B.1.1.7),
                               diffExpLFCshrink_vero24h %>% 
                                 dplyr::select(gene,log2FoldChange) %>% 
                                 dplyr::rename("SARS-CoV-2 infected Vero 24h"=log2FoldChange),
                               by="gene")

lfc_B.1.351vsVero24h <- inner_join(log2FC.all.df %>% dplyr::select(gene,B.1.351),
                                   diffExpLFCshrink_vero24h %>% 
                                     dplyr::select(gene,log2FoldChange) %>% 
                                     dplyr::rename("SARS-CoV-2 infected Vero 24h"=log2FoldChange),
                                   by="gene")

lfc_SamplesvsVero24h_vero8h <- inner_join(log2FC.all.df,
                                   diffExpLFCshrink_vero24h %>% 
                                     dplyr::select(gene,log2FoldChange) %>% 
                                     dplyr::rename("SARS-CoV-2 infected Vero 24h"=log2FoldChange),
                                   by="gene")

lfc_SamplesvsVero24h_vero8h <- inner_join(lfc_SamplesvsVero24h_vero8h,
                                      diffExpLFCshrink_vero8h %>% 
                                        dplyr::select(gene,log2FoldChange) %>% 
                                        dplyr::rename("SARS-CoV-2 infected Vero 8h"=log2FoldChange),
                                      by="gene")

lfc_SamplesvsVero24h_vero8h <- lfc_SamplesvsVero24h_vero8h %>% filter(!is.na(`SARS-CoV-2 infected Vero 8h`)) %>%  
  filter(!is.na(`SARS-CoV-2 infected Vero 24h`)) %>% 
  dplyr::select(-"SARS-CoV-2 infected Vero 8h")

limitRange <- function(data, mapping, ...) { 
  ggplot(data = data, mapping = mapping, ...) + 
    geom_point(...,size=0.9,alpha=0.3) + 
    geom_smooth(method = "lm", se = FALSE) +
    scale_y_continuous(limits = c(-6, 6)) +
    scale_x_continuous(limits = c(-6, 6)) 
}


ggpairs(lfc_SamplesvsVero24h_vero8h, columns=2:5,diag ="blankDiag", 
        upper=list(continuous=wrap("cor",method="pearson", size=5)), 
        #lower=list(continuous=wrap("points",alpha=0.3,size=0.9)), 
        lower=list(continuous=limitRange), 
        columnSamplesels = str_wrap(colnames(lfc_SamplesvsVero24h_vero8h)[-1],20)) -> pairedCorPlotPublic

pairedCorPlotPublic + 
  labs(x=str_wrap("log2 Fold Change \n (SARS-CoV-2 Infected / Non-infected)",30),
       y=str_wrap("log2 Fold Change \n (SARS-CoV-2 Infected / Non-infected)",30)) + 
  theme_bw(base_size = 10) + theme(strip.text = element_text(size=10)) -> pairedCorPlotPublic

ggsave("/figures/correlationWithPublicData_labSamples_log2FoldChanges_matrix.pdf", width = 12,height = 12, plot = pairedCorPlotPublic)


svg("/figures/correlationSamples_B.1_B.1.1.7_B.1.351_log2FoldChanges.svg", width = 6,height = 6)

corrplot(cor(lfc_SamplesvsVero24h_vero8h[,c(2:4)],),
         addrect = 0, 
         addCoef.col = 'white',
         method = 'circle',
         type="upper",
         number.cex = 0.75,tl.cex = 1,
         tl.srt = 45,
         #         col = rev(coolwarm(100)),
         diag = T,
         is.corr = T,tl.col = "black")

dev.off()


signif24h_topUp <- diffExpLFCshrink_vero24h %>% filter(padj < 0.1 &! is.na(padj) & !(gene %in% cov2genes) &! grepl("gene-",gene)) %>% arrange(desc(log2FoldChange)) %>% dplyr::select(gene) %>% distinct() %>% unlist()

signif24h_topDown <- diffExpLFCshrink_vero24h %>% filter(padj < 0.1 &! is.na(padj) & !(gene %in% cov2genes) &! grepl("gene-",gene)) %>% arrange(log2FoldChange) %>% dplyr::select(gene) %>% distinct() %>% unlist()
#%>% head(200)


commondiffexpGenesAllStrains <- intersect(intersect(diffexpGenes_wu,diffexpGenes_uk),diffexpGenes_sa)
  #intersect(intersect(intersect(diffexpGenes_wu,diffexpGenes_uk),diffexpGenes_sa),c(signif24h_topDown,signif24h_topUp))
#uniondiffexpGenesAllStrains <- union(union(diffexpGenes_wu,diffexpGenes_uk),diffexpGenes_sa)
ggpairs(lfc_SamplesvsVero24h_vero8h %>% 
          #filter(gene %in% commondiffexpGenesAllStrains), 
          #filter(gene %in% uniondiffexpGenesAllStrains), 
        filter(gene %in% c(signif24h_topDown,signif24h_topUp)), 
        columns = c(2:5),
        diag ="blankDiag", 
        upper=list(continuous=wrap("cor",method="pearson", size=5)), 
        lower=list(continuous=limitRange)
        )
 #%>% filter(gene %in% diffexpGenes_wu)
#%>% filter(gene %in% c(signif24h_topDown,signif24h_topUp))
#aes(label=paste(..rr.label.., ..p.label..,sep="~`,`~")

ngenes <- 300
  top50up_wu <- res_all_lfcShrink %>% filter(padj < 0.01 &! is.na(padj) & contrast=="wu_vs_vero" & log2FoldChange > 1 & !(gene %in% cov2genes) & !grepl("(gene|rna)",gene,perl = T)) %>% dplyr::arrange(desc(log2FoldChange)) %>% head(ngenes) %>% dplyr::select(gene) %>% distinct() %>% unlist()
  
  top50dw_wu <- res_all_lfcShrink %>% filter(padj < 0.01 &! is.na(padj) & contrast=="wu_vs_vero" & log2FoldChange < -1 & !(gene %in% cov2genes) &  !grepl("(gene|rna)",gene,perl = T)) %>% dplyr::arrange(log2FoldChange)  %>% head(ngenes) %>% dplyr::select(gene) %>% distinct() %>% unlist()
  
  top50up_uk <- res_all_lfcShrink %>% filter(padj < 0.01 &! is.na(padj) & contrast=="uk_vs_vero" & log2FoldChange > 1 & !(gene %in% cov2genes) & !grepl("(gene|rna)",gene,perl = T)) %>% dplyr::arrange(desc(log2FoldChange)) %>% head(ngenes) %>% dplyr::select(gene) %>% distinct() %>% unlist()
  
  top50dw_uk <- res_all_lfcShrink %>% filter(padj < 0.01 &! is.na(padj) & contrast=="uk_vs_vero" & log2FoldChange < -1 & !(gene %in% cov2genes) & !grepl("(gene|rna)",gene,perl = T)) %>% dplyr::arrange(log2FoldChange) %>% head(ngenes) %>% dplyr::select(gene) %>% distinct() %>% unlist()
  
  top50up_sa <- res_all_lfcShrink %>% filter(padj < 0.01 &! is.na(padj) & contrast=="sa_vs_vero" & log2FoldChange > 1 & !(gene %in% cov2genes) & !grepl("(gene|rna)",gene,perl = T)) %>% dplyr::arrange(desc(log2FoldChange)) %>% head(ngenes) %>% dplyr::select(gene) %>% distinct() %>% unlist()
  
  top50dw_sa <- res_all_lfcShrink %>% filter(padj < 0.01 &! is.na(padj) & contrast=="sa_vs_vero" & log2FoldChange < -1 & !(gene %in% cov2genes) & !grepl("(gene|rna)",gene,perl = T)) %>% dplyr::arrange(log2FoldChange) %>% head(ngenes) %>% dplyr::select(gene) %>% distinct() %>% unlist()


lfc_SamplesvsVero24h_vero8h %>%
  dplyr::select(gene,B.1,"SARS-CoV-2 infected Vero 24h") %>% 
  filter(gene %in% diffexpGenes_wu) %>% 
  ggplot(aes(`SARS-CoV-2 infected Vero 24h`,B.1)) + 
  geom_point(size=1,alpha=0.5) + 
  geom_smooth(method = "lm", se = TRUE) + 
  stat_cor(method = "pearson",label.x = -5,label.y = 5,size=3.5) + 
  coord_cartesian(xlim=c(-6,6),ylim = c(-6,6)) +
  theme(aspect.ratio = 0.5) +
  theme_prism(base_size = 10) +
  labs(y="log2 Fold Change \n (Infected / Non-Infected) \n B.1",
       x="log2 Fold Change \n (Infected / Non-Infected) \n Riva et al. 2020") -> pCor24_wu

lfc_SamplesvsVero24h_vero8h %>% 
  dplyr::select(gene,B.1.1.7,"SARS-CoV-2 infected Vero 24h") %>% 
  filter(gene %in% diffexpGenes_uk) %>%
  ggplot(aes(`SARS-CoV-2 infected Vero 24h`,B.1.1.7)) + 
  geom_point(size=1,alpha=0.5) + 
  geom_smooth(method = "lm", se = TRUE) + 
  stat_cor(method = "pearson",label.x = -5,label.y = 5,size=3.5) + 
  coord_cartesian(xlim=c(-6,6),ylim = c(-6,6)) + 
  theme(aspect.ratio = 0.5) +
  theme_prism(base_size = 10) +
  labs(y="log2 Fold Change \n (Infected / Non-Infected) \n B.1.1.7",
       x="log2 Fold Change \n (Infected / Non-Infected) \n Riva et al. 2020") -> pCor24_uk

lfc_SamplesvsVero24h_vero8h %>% 
  dplyr::select(gene,B.1.351,"SARS-CoV-2 infected Vero 24h") %>% 
  filter(gene %in% diffexpGenes_sa) %>%
  ggplot(aes(`SARS-CoV-2 infected Vero 24h`,B.1.351)) + 
  geom_point(size=1,alpha=0.5) + 
  geom_smooth(method = "lm", se = TRUE) + 
  stat_cor(method = "pearson",label.x = -5,label.y = 5, size=3.5) + 
  coord_cartesian(xlim=c(-6,6),ylim = c(-6,6)) + 
  theme(aspect.ratio = 0.5) +
  theme_prism(base_size = 10) +
labs(y="log2 Fold Change \n (Infected / Non-Infected) \n B.1.351",
     x="log2 Fold Change \n (Infected / Non-Infected) \n Riva et al. 2020") -> pCor24_sa

  p <- pCor24_wu + pCor24_uk + pCor24_sa 
  p
ggsave("/figures/correlationPublicData_vero24hRiva_diffExpGenes_byStrain.svg",width = 14, height = 4, plot = p)

#%>% filter(gene %in% commondiffexpGenesAllStrains)

#range01 <- function(x){(x-min(x))/(max(x)-min(x))}

# lfc_SamplesvsVero24h_vero8h %>% mutate(z_B.1=(B.1-median(B.1))/sd(B.1),
#                                    z_B.1.1.7=(B.1.1.7-median(B.1.1.7))/sd(B.1.1.7),
#                                    z_B.1.351=(B.1.351-median(B.1.351))/sd(B.1.351),
#                                    z_24h=(`SARS-CoV-2 infected Vero 24h`-median(`SARS-CoV-2 infected Vero 24h`))/sd(`SARS-CoV-2 infected Vero 24h`))
#%>% summarise(range(z_B.1),range(z_B.1.1.7),range(z_B.1.351),range(z_24h))

#unique(c(top50up_wu,top50dw_wu,top50up_uk,top50dw_uk,top50up_sa,top50dw_sa))
lfcToPubData.m <- lfc_SamplesvsVero24h_vero8h %>% 
  filter(gene %in% c(intersect(intersect(top50dw_sa,top50dw_uk),top50dw_wu)[1:30],
                     intersect(intersect(top50up_sa,top50up_uk),top50up_wu)[1:30])) %>% 
  dplyr::select(-gene) %>% as.matrix()

# %>% 
#   mutate(z_B.1=(B.1-median(B.1))/sd(B.1),
#          z_B.1.1.7=(B.1.1.7-median(B.1.1.7))/sd(B.1.1.7),
#          z_B.1.351=(B.1.351-median(B.1.351))/sd(B.1.351),
#          z_24h=(`SARS-CoV-2 infected Vero 24h`-median(`SARS-CoV-2 infected Vero 24h`))/sd(`SARS-CoV-2 infected Vero 24h`)) %>% 
#   filter(gene %in% unique(c(top50up_wu,top50dw_wu,top50up_uk,top50dw_uk,top50up_sa,top50dw_sa))) %>% dplyr::select(dplyr::starts_with("z"),-gene) %>% as.matrix()
#%>% filter(gene %in% commondiffexpGenesAllStrains)
rownames(lfcToPubData.m) <- lfc_SamplesvsVero24h_vero8h %>% 
  filter(gene %in% c(intersect(intersect(top50dw_sa,top50dw_uk),top50dw_wu)[1:30],
                     intersect(intersect(top50up_sa,top50up_uk),top50up_wu)[1:30])) %>% 
  dplyr::select(gene) %>% unlist()
  #filter(gene %in% unique(c(top50up_wu,top50dw_wu,top50up_uk,top50dw_uk,top50up_sa,top50dw_sa))) 

colnames(lfcToPubData.m) <- str_wrap(colnames(lfcToPubData.m),20)
myscale <- c(rev(brewer.blues(51)[10:50]),"white",brewer.reds(51)[10:50])
mybreaks <- c(seq(-6,-2,by=0.1),0,seq(2,6,by=0.1))
pheatmap(t(lfcToPubData.m), color = coolwarm(100), border_color = NA,
         angle_col = 45,fontsize_col = 6,cellwidth = 10,cellheight = 25,treeheight_row = 25,treeheight_col = 15,)

#ggpairs(lfc_SamplesvsVero24h_vero8h %>% filter(gene %in% uniondiffexpGenesAllStrains), columns = c(2:5))


#patientLFC <- openxlsx::read.xlsx("/plots/cov2_figures_2022/differentially_expressed_genes_patients_covid19_Blanco-Melo2020.xlsx",sheet = 2)

#ggpairs(inner_join(lfc_SamplesvsVero24h_vero8h,patientLFC %>% filter(log2FoldChange!= "NaN") %>% dplyr::select(Gene_name,log2FoldChange) %>% rename("Patients_Melo"=log2FoldChange) %>% mutate(Patients_Melo=as.numeric(Patients_Melo)),by=c("gene"="Gene_name")), columns = 2:6)

#----- comparison to GSE30589, fetched using GEO2R script
#https://www.ncbi.nlm.nih.gov/geo/geo2r/?acc=GSE30589
#ex.df %>% openxlsx::write.xlsx(., "/tables/geo2R_GSE30589_expression.xlsx")
ex.df %>% 


#--------------------------------Correlation Matrix log2 Fold Changes ------------------------------·
lfc_SamplesvsVero24h_vero8h.m <- as.matrix(lfc_SamplesvsVero24h_vero8h[,-1])
#cor_plot(cor(lfc_SamplesvsVero24h_vero8h.m),label = T)

svg("/figures/correlationWithPublicData_labSamples_log2FoldChanges.svg", width = 12,height = 12)

corrplot(cor(lfc_SamplesvsVero24h_vero8h.m),
order = 'hclust', 
         addrect = 0, 
         addCoef.col = 'white',
         method = 'circle',
         type="upper",
         number.cex = 0.6,
         tl.srt = 45,
         #         col = rev(coolwarm(100)),
         diag = T,
         is.corr = T,tl.col = "black")

#colorlegend(xlim=c(10,15), ylim=c(10,15), rev(coolwarm(100)), c(seq(-1,1,.25)), align="l", vertical=TRUE, addlabels=TRUE)

dev.off()



#------------------------------------

cor_lfc_B.1vsVero8h <- lfc_B.1vsVero8h %>% filter(gene %in% signif_vero8h) %>% 
  cor_test(method = "pearson", vars = c("B.1","SARS-CoV-2 infected Vero 8h")) %>% adjust_pvalue()

lfc_B.1vsVero8h %>% filter(gene %in% signif_vero8h) %>% pairedCor2(.,labelx = "B.1",labely = "SARS-CoV-2 infected \n Zhang Y et al. 2021") 
#+ annotate(geom = "text", x=-2, y=3,label=paste0("R",cor_lfc_B.1vsVero8h$cor)

inner_join(log2FC.all.df %>% dplyr::select(gene,B.1.1.7),
           diffExpLFCshrink_vero8h %>% 
  dplyr::select(gene,log2FoldChange) %>% 
  dplyr::rename("SARS-CoV-2 infected Vero 8h"=log2FoldChange),
by="gene") %>% pairedCor(.,labelx = "B.1.1.7",labely = "SARS-CoV-2 infected Vero 8h")

inner_join(log2FC.all.df %>% dplyr::select(gene,B.1.351),
           diffExpLFCshrink_vero8h %>% 
             dplyr::select(gene,log2FoldChange) %>% 
             dplyr::rename("SARS-CoV-2 infected Vero 8h"=log2FoldChange),
           by="gene") %>% pairedCor(.,labelx = "B.1.351",labely = "SARS-CoV-2 infected Vero 8h")




inner_join(log2FC.all.df %>% dplyr::select(gene,B.1),
           diffExpLFCshrink_vero24h %>% 
             dplyr::select(gene,log2FoldChange) %>% 
             dplyr::rename("SARS-CoV-2 infected Vero 24h"=log2FoldChange),
           by="gene") %>% pairedCor(.,labelx = "B.1",labely = "SARS-CoV-2 infected \n Riva L et al. 2020")

inner_join(log2FC.all.df %>% dplyr::select(gene,B.1.1.7),
           diffExpLFCshrink_vero24h %>% 
             dplyr::select(gene,log2FoldChange) %>% 
             dplyr::rename("SARS-CoV-2 infected Vero 24h"=log2FoldChange),
           by="gene") %>% pairedCor(.,labelx = "B.1.1.7",labely = "SARS-CoV-2 infected \n Riva L et al. 2020")

inner_join(log2FC.all.df %>% 
             filter(gene %in% signif_vero8h) %>% 
             dplyr::select(gene,B.1.351),
           diffExpLFCshrink_vero24h %>% 
             filter(gene %in% signif_vero8h) %>%
             dplyr::select(gene,log2FoldChange) %>% 
             dplyr::rename("SARS-CoV-2 infected Vero 24h"=log2FoldChange),
           by="gene") %>% pairedCor(.,labelx = "B.1.351",labely = "SARS-CoV-2 infected \n Riva L et al. 2020")






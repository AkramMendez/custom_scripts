#!/usr/bin/env Rscript

library("Rsubread")
library("DESeq2")
library("tidyverse")
library("pheatmap")
library("RColorBrewer")

threads=1
inputdir_rev="/mappings/infected_cells_mapped_to_sarscov2/infected_nodup_cov2only/rev/"

bamfiles_rev<-list.files(inputdir_rev,pattern = ".bam$")
bamfiles_rev<-paste0(inputdir_rev,bamfiles_rev)
bamfiles_rev<-bamfiles_rev[!grepl("Vero",bamfiles_rev)]



#annotation_file<-"chlSab2_wuhCor1_ecoliK12_ncbiGenes_chroms_concat.gtf"
annotation_file="~/genomes/sarsCov2/genes/ncbiGenes.gtf"
metadata_rev<- NULL
metadata_rev$sample<-as.factor(basename(bamfiles_rev))
metadata_rev$condition<-as.factor(c(rep("infected",9),rep("uninfected",0)))
metadata_rev<-as.data.frame(metadata_rev)
rownames(metadata_rev)<-basename(bamfiles_rev)



featureCounts_data_rev<-featureCounts(bamfiles_rev, annot.ext=annotation_file, isGTFAnnotationFile = TRUE, nthreads = threads)

all(colnames(featureCounts_data_rev$counts) %in% rownames(metadata_rev))
dds_rev<-DESeqDataSetFromMatrix(countData = featureCounts_data_rev$counts, colData=metadata_rev, design= ~ sample)

dds_rev<- estimateSizeFactors(dds_rev)

sizeFactors(dds_rev)


normalized_counts_rev <- counts(dds_rev, normalized=F)

colnames(normalized_counts_rev)<-gsub("(_S\\d+.*.sorted.|sarscov2only.bam|renamed.hisat2.nodup.uniq.sorted.sarscov2Only.rev.bam)","",colnames(normalized_counts_rev),perl = T)

normalized_counts_rev %>% as.data.frame() %>% mutate(gene=rownames(normalized_counts_rev)) %>% gather(.,key="sample",value = "norm.counts",-gene) %>% group_by(sample) %>% summarise(total_norm_counts=sum(sum(norm.counts))) %>% as.data.frame() -> norm.counts.rev.df

norm.counts.rev.df %>% ggplot(.,aes(sample,total_norm_counts)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle=60, hjust = 1)) + labs(title = "Norm Counts Rev")


inputdir_fwd="/mappings/infected_cells_mapped_to_sarscov2/infected_nodup_cov2only/fwd/"

bamfiles_fwd<-list.files(inputdir_fwd,pattern = ".bam$")
bamfiles_fwd<-paste0(inputdir_fwd,bamfiles_fwd)
bamfiles_fwd<-bamfiles_fwd[!grepl("Vero",bamfiles_fwd)]

metadata_fwd<- NULL
metadata_fwd$sample<-as.factor(basename(bamfiles_fwd))
metadata_fwd$condition<-as.factor(c(rep("infected",9),rep("uninfected",0)))
metadata_fwd<-as.data.frame(metadata_fwd)
rownames(metadata_fwd)<-basename(bamfiles_fwd)



featureCounts_data_fwd<-featureCounts(bamfiles_fwd, annot.ext=annotation_file, isGTFAnnotationFile = TRUE)

all(colnames(featureCounts_data_fwd$counts) %in% rownames(metadata_fwd))
dds_fwd<-DESeqDataSetFromMatrix(countData = featureCounts_data_fwd$counts, colData=metadata_fwd, design= ~ sample)

dds_fwd<- estimateSizeFactors(dds_fwd)

sizeFactors(dds_fwd)


normalized_counts_fwd <- counts(dds_fwd, normalized=F)

colnames(normalized_counts_fwd)<-gsub("(_S\\d+.*.sorted.|sarscov2only.bam|renamed.hisat2.nodup.uniq.sorted.sarscov2Only.fwd.bam)","",colnames(normalized_counts_fwd),perl = T)

normalized_counts_fwd %>% as.data.frame() %>% mutate(gene=rownames(normalized_counts_fwd)) %>% gather(.,key="sample",value = "norm.counts",-gene) %>% group_by(sample) %>% summarise(total_norm_counts=sum(sum(norm.counts))) %>% as.data.frame() -> norm.counts.fwd.df

norm.counts.fwd.df %>% ggplot(.,aes(sample,total_norm_counts)) + geom_bar(stat="identity") + theme(axis.text.x = element_text(angle=60, hjust = 1)) + labs(title = "Norm Counts Fwd")

#__________________________________________________________#

# No normalized counts:
colnames(featureCounts_data_fwd$counts)<-gsub("(_S\\d+.*.sorted.|sarscov2only.bam|renamed.hisat2.nodup.uniq.sorted.sarscov2Only.|.sarscov2Only.fwd.bam|_renamed.hisat2.nodup.uniq.sorted.sarscov2Only.fwd.bam)","",colnames(featureCounts_data_fwd$counts),perl = T)
colnames(featureCounts_data_rev$counts)<-gsub("(_S\\d+.*.sorted.|sarscov2only.bam|renamed.hisat2.nodup.uniq.sorted.sarscov2Only.|.sarscov2Only.rev.bam|_renamed.hisat2.nodup.uniq.sorted.sarscov2Only.rev.bam)","",colnames(featureCounts_data_rev$counts),perl = T)

featureCounts_data_fwd$counts %>% as.data.frame() %>% mutate(gene=rownames(featureCounts_data_fwd$counts)) %>% gather(.,key="sample",value = "raw_counts",-gene) %>% group_by(sample) %>% summarise(total_raw_counts=sum(sum(raw_counts))) %>% as.data.frame() -> total_raw_counts_fwd
featureCounts_data_rev$counts %>% as.data.frame() %>% mutate(gene=rownames(featureCounts_data_rev$counts)) %>% gather(.,key="sample",value = "raw_counts",-gene) %>% group_by(sample) %>% summarise(total_raw_counts=sum(sum(raw_counts))) %>% as.data.frame() -> total_raw_counts_rev
colors<-c("#66c2a4","#2ca25f","#006d2c","#6baed6","#3182bd","#08519c","#fb6a4a","#de2d26","#a50f15")


total_raw_counts_fwd %>% 
  mutate(sample=gsub("_Infected_cell_"," Infected ",sample)) %>%
  mutate(sample=if_else(
    grepl("COVID_19",sample),gsub("COVID_19","WU",sample),
    if_else(grepl("COVID_SA",sample),gsub("COVID_SA","SA",sample),
            if_else(grepl("COVID_UK",sample),gsub("COVID_UK","UK",sample),"NA")))) %>% 
  ggplot(.,aes(sample,total_raw_counts,fill=sample)) + 
  geom_bar(stat = "identity", width = 0.5) + 
  scale_fill_manual(values = colors) + theme_bw() + theme(axis.text.x = element_text(angle=60,hjust = 1)) + xlab("") + ylab("Total counts") + labs(title = "Total counts",subtitle= "SARS-Cov2 Forward",fill="") ->plot_total_fwd
ggplotly(plot_total_fwd)

total_raw_counts_rev %>% mutate(sample=gsub("_Infected_cell_"," Infected ",sample)) %>%
  mutate(sample=if_else(
    grepl("COVID_19",sample),gsub("COVID_19","WU",sample),
    if_else(grepl("COVID_SA",sample),gsub("COVID_SA","SA",sample),
            if_else(grepl("COVID_UK",sample),gsub("COVID_UK","UK",sample),"NA")))) %>% 
  ggplot(.,aes(sample,total_raw_counts,fill=sample)) + 
  geom_bar(stat = "identity", width = 0.5) + 
  scale_fill_manual(values = colors) + theme_bw() + theme(axis.text.x = element_text(angle=60,hjust = 1)) + xlab("") + ylab("Total counts") + labs(title = "Total counts",subtitle= "SARS-Cov2 Forward",fill="") ->plot_total_rev
ggplotly(total_raw_counts_rev)


#Raw count plots by Gene
left_join(featureCounts_data_fwd$counts %>% as.data.frame() %>% mutate(gene=rownames(featureCounts_data_fwd$counts)) %>% gather(.,key="sample",value = "raw_counts",-gene),total_raw_counts_fwd,by=c("sample")) %>% mutate(percentage_raw_counts=(raw_counts/total_raw_counts)*100) %>% ggplot(.,aes(gene,raw_counts, color=sample)) + geom_point(size=1,alpha=0.9) + coord_flip() + theme_bw() + scale_color_manual(values = colors) + ylab("Counts") + xlab("Gene") + labs(title="Counts by gene Fwd",subtitle = "SARS-Cov2 Forward",color="") ->plot_fwd

left_join(featureCounts_data_rev$counts %>% as.data.frame() %>% mutate(gene=rownames(featureCounts_data_rev$counts)) %>% gather(.,key="sample",value = "raw_counts",-gene),total_raw_counts_rev,by=c("sample")) %>% mutate(percentage_raw_counts=(raw_counts/total_raw_counts)*100) %>% ggplot(.,aes(gene,raw_counts, color=sample)) + geom_point(size=1,alpha=0.9) + coord_flip() + theme_bw() + scale_color_manual(values = colors) + ylab("Counts") + xlab("Gene") + labs(title="Counts by gene Rev",subtitle = "SARS-Cov2 Reverse",color="") ->plot_rev

ggplotly(plot_fwd)

ggplotly(plot_rev)
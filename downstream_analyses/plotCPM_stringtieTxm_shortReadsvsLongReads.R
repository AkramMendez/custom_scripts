library(rtracklayer)
library(tidyverse)
library(data.table)


colNames <- c("chr","source","type","start","end","score","strand","frame","attribute","TPM")

cpmTxm_shortReads <- fread("bedtools_intersect_extractedTPM_stringtie2_shortReadsTxm_toT2Tannotation.tsv", sep="\t", col.names = c(paste0(colNames,"_shortTxm"),paste0(colNames[-10],"_TxmT2T"))) %>% filter(type_TxmT2T=="transcript") %>% separate(attribute_TxmT2T,into = c("transcript_id","gene_id","gene_name"), sep = ";") %>% mutate(transcript_id=gsub("\\s{0,}transcript_id\\s+","",transcript_id),gene_id=gsub("\\s{0,}gene_id\\s+","",gene_id), gene_name=gsub("\\s{0,}gene_name\\s+","",gene_name))

cpmTxm_longReads <- fread("Control_bedtools_intersect_extractedTPM_stringtie2_longReadsTxm_toT2Tannotation.tsv", sep="\t", col.names = c(paste0(colNames,"_longTxm"),paste0(colNames[-10],"_TxmT2T"))) %>% filter(type_TxmT2T=="transcript") %>% separate(attribute_TxmT2T,into = c("transcript_id","gene_id","gene_name"), sep = ";") %>% mutate(transcript_id=gsub("\\s{0,}transcript_id\\s+","",transcript_id),gene_id=gsub("\\s{0,}gene_id\\s+","",gene_id), gene_name=gsub("\\s{0,}gene_name\\s+","",gene_name))



group_by(GENEID) %>% top_n(1,wt = width) %>% filter(., rank(dplyr::desc(width), ties.method = "first")==1)


cpmTxm_shortReads_df <- cpmTxm_shortReads %>% dplyr::select(attribute_shortTxm,TPM_shortTxm,gene_name) %>% distinct() %>%  group_by(gene_name) %>% mutate(meanTPM=mean(TPM_shortTxm)) %>% ungroup() %>% dplyr::select(gene_name,meanTPM) %>% distinct() %>% mutate(gene_name=gsub("\"","",gene_name))


cpmTxm_longReads_df <- cpmTxm_longReads %>% dplyr::select(attribute_longTxm,TPM_longTxm,gene_name) %>% distinct() %>%  group_by(gene_name) %>% mutate(meanTPM=mean(TPM_longTxm)) %>% ungroup() %>% dplyr::select(gene_name,meanTPM) %>% distinct() %>% mutate(gene_name=gsub("\"","",gene_name))


df <- inner_join(cpmTxm_shortReads_df,cpmTxm_longReads_df, by="gene_name") %>% replace(is.na(.),0) %>% rename(meanTPM_shortReads=meanTPM.x,meanTPM_longReads=meanTPM.y)


df %>% ggplot(aes(log10(meanTPM_shortReads + 1),log10(meanTPM_longReads + 1))) + geom_point() + stat_regline_equation(aes(label=..rr.label..)) + theme_prism()



# featureCounts -----------------------------------------------------------

libsize_control_ont <- 2325526
libsize_control_u2os <- 140851358


calculateTPM <- function(df, libSize){
  
  res <- df %>% 
    mutate(rpk=counts/(length/1000)) %>%
    mutate(tpm=(rpk/sum(rpk))*10^6) %>% 
    ungroup()
}

colsFeatCounts <-c("gene_name","chr","start","end","strand","length","counts")

countsT2T_short <- fread("featureCounts_Control_input_spikeInNorm_ILMN_mappedToGenomeT2T.tsv", col.names = colsFeatCounts,sep="\t") %>% mutate(sample="Control_short")

countsT2T_long <- fread("featureCounts_Control_ONT_mappedToGenomeT2T.tsv", col.names = colsFeatCounts,sep="\t") %>% mutate(sample="Control_long")


countsT2T_short <- countsT2T_short %>% 
mutate(rpk=counts/(length/1000)) %>% #Calculate reads per kilobase
  mutate(tpm=(rpk/sum(rpk))*10^6) #Calculate TPM


countsT2T_long <- countsT2T_long %>% 
  mutate(rpk=counts/(length/1000)) %>% #Calculate reads per kilobase
  mutate(tpm=(rpk/sum(rpk))*10^6) #Calculate TPM


df <- inner_join(countsT2T_short %>% dplyr::select(gene_name,counts,tpm_short=tpm,sample),
            countsT2T_long %>% dplyr::select(gene_name,counts,tpm_long=tpm,sample),
            by="gene_name"
)

#filter(tpm.x >= 0,tpm.y >= 0)
df %>% mutate(x=log10(tpm_short + 1), y=log10(tpm_long + 1)) %>% 
  ggplot(aes(x,y)) + geom_point(size=0.01, alpha=0.5) + geom_smooth(method = "lm", se = F) + 
  stat_cor(method = "spearman",na.rm = T) + theme_prism() + 
  scale_y_continuous(guide = "prism_offset_minor", expand = expansion(mult = c(0.01,NA))) + 
  scale_x_continuous(guide = "prism_offset_minor", expand = expansion(mult = c(0.01,NA))) +
  theme(aspect.ratio = 1) +
  labs(title = "Correlation expression \n Nanopore long-reads vs Illumina short-reads U2OS \n TPM", 
       x="log10(TPM + 1) \n Expression short-reads",
       y="Expression long-reads \n log10(TPM + 1)") -> p


ggsave("figure_compilation/correlation_TPMexpression_Nanopore_vs_Illumina_U2OS.pdf", width = 7, height = 7)


######### Quantification using Salmon (for short reads) and Nanocount (long-reads)


nanocount_txm <- fread("/nanoCount/nanoCount_ONT_minimap2.24_t2tTxm.chm13.txm.14.mapped.tsv", header = T,sep="\t")  %>% mutate(sample="Long_reads")

salmon_control <- fread("quant.sf",header = T,sep="\t") %>% mutate(sample="Short_reads")

df <- inner_join(salmon_control %>% dplyr::select(Name,tpm_short=TPM,sample),
                 nanocount_txm %>% dplyr::select(transcript_name,tpm_long=tpm,sample),
                 by=c("Name"="transcript_name")
)



df %>% as.data.frame() %>% mutate(x=log10(tpm_short + 0.01), y=log10(tpm_long  + 0.01)) %>% 
  ggplot(aes(x,y)) + geom_point(size=0.01, alpha=0.5) + geom_smooth(method = "lm", se = F) + 
  stat_cor(method = "spearman",na.rm = T) + theme_prism() + 
  scale_y_continuous(guide = "prism_offset_minor", expand = expansion(mult = c(0.01,NA))) + 
  scale_x_continuous(guide = "prism_offset_minor", expand = expansion(mult = c(0.01,NA))) +
  theme(aspect.ratio = 1) +
  labs(title = "Correlation expression \n Nanopore long-reads vs Illumina short-reads \n TPM", 
       x="log10(TPM + 0.01) \n Expression short-reads",
       y="Expression long-reads \n log10(TPM + 0.01)") #-> p

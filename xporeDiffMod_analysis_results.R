library(tidyverse)
library(data.table)
library(plyranges)
library(rtracklayer)


chrEndCoords <- rtracklayer::import("~/gitrepos/bioinfopipes/lab/CHM13_hg38chrY_chromEndCoordinates_50kb.genome",format="BED")

chrEndCoords <- chrEndCoords %>% as.data.frame() %>% mutate(strand=if_else(end > 50000,"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)

xptable <- fread("~/MondalLab/terra/nanoporeTerra/xpore_diffmod/diffmod.table", sep = ",",header = T)

gff <- rtracklayer::import("/home/akram/MondalLab/terra/CHM13/CHM13_hg38chrY_ecoliK12_gffreadFormatted.genePred")

gff <- gff %>% filter(seqnames!="U00096.2")

gff <- gff %>% as.data.frame() %>% dplyr::select(seqnames,start,end,width,strand,gene_id,gene_name) %>% distinct()
genePred <- fread("~/MondalLab/terra/CHM13/CHM13_hg38chrY_ecoliK12_gffreadFormatted.genePred",sep = "\t", header = F)[,c(1:5,12)]

colnames(genePred) <- c("tx_id","chr","gene_strand","start","end","gene_id")
genePred <- genePred %>% filter(chr!="U00096.2")

genePred <- genePred %>% 
  group_by(gene_id) %>% 
  mutate(gene_start=min(start),gene_end=max(end)) %>%
  dplyr::select(chr,gene_start,gene_end,gene_strand,gene_id) %>%
  distinct() %>% 
  ungroup()

genePred <- genePred %>% mutate(gene_width=gene_end-gene_start)

xptable <- left_join(xptable,genePred,by=c("id"="gene_id"))
#"(A|G|T)(A|G)AC(A|C|T)"
xptable <- xptable %>% mutate(drach= if_else(grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T),TRUE,FALSE))
xptable <- xptable %>% filter(!is.na(gene_start))
xptable <- xptable %>% mutate(start=position-3,end=position+2) %>% mutate(width=end-start)

#%>% dplyr::select(chr,gene_start,gene_end,gene_width,gene_strand,colnames(xptable)[1:19],drach)
xptable.gr <- xptable %>% dplyr::select(chr,start,end,width,gene_strand,colnames(xptable)[1:19],drach) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T,seqnames.field = "chr",strand.field = "gene_strand",end.field = "end",start.field = "start")

join_overlap_inner_directed(xptable.gr,chrEndCoords)


xptable.gr %>% as.data.frame() %>% dplyr::select(kmer,mod_rate_CTRLSH.REP1,position) %>% filter(grepl("((..AC.)|GGATT|GAAGG|TAAAG|TAAGA)",kmer,perl = T)) %>% arrange(dply::desc(mod_rate_CTRLSH.REP1)) %>% head(100)


xptable %>% filter(diff_mod_rate_CTRLSH_vs_METTL3KD > 0 & pval_CTRLSH_vs_METTL3KD < 0.05 & grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T) ) %>% dplyr::select(kmer,id,position,starts_with("mod_rate")) %>% pivot_longer(names_to = "cond",values_to = "mod_rate", -c("kmer","id","position")) %>% ggplot(aes(mod_rate,fill=cond)) + geom_density(alpha=0.5)


#(G|A|T)GT(T|C)(T|C)
modRatesAtChrEnds_csh <- join_overlap_inner_directed(xptable.gr,chrEndCoords)

modRatesAtChrEnds_csh

#modRatesAtChrEnds_csh 
xptable.gr %>% as.data.frame() %>% filter(mod_rate_CTRLSH.REP1 > 0 & grepl("..A..",kmer,perl = T)) %>% dplyr::select(seqnames,start,end,kmer,mod_rate_CTRLSH.REP1,strand) %>% fwrite(.,"~/MondalLab/terra/nanoporeTerra/xpore_modRate_NNANN_ctrlsh.bed", sep="\t",col.names = F)

xptable.gr %>% as.data.frame() %>% filter(mod_rate_CTRLSH.REP1 > 0 & grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T)) %>% dplyr::select(seqnames,start,end,kmer,mod_rate_CTRLSH.REP1,strand) %>% fwrite(.,"~/MondalLab/terra/nanoporeTerra/xpore_modRate_DRACH_ctrlsh.bed", sep="\t",col.names = F)

modRatesAtChrEnds_mettl3kd <- join_overlap_inner_directed(xptable.gr,chrEndCoords)


#modRatesAtChrEnds_mettl3kd 

xptable.gr %>% as.data.frame() %>% filter(mod_rate_METTL3KD.REP1 > 0 & grepl("..A..",kmer,perl = T)) %>% mutate(width=end-start) %>% dplyr::select(seqnames,start,end,kmer,mod_rate_METTL3KD.REP1,strand) %>% fwrite(.,"~/MondalLab/terra/nanoporeTerra/xpore_modRate_NNANN_mettl3kd.bed", sep="\t",col.names = F)

xptable.gr %>% as.data.frame() %>% filter(mod_rate_METTL3KD.REP1 > 0 & grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T)) %>% mutate(width=end-start) %>% dplyr::select(seqnames,start,end,kmer,mod_rate_METTL3KD.REP1,strand) %>% fwrite(.,"~/MondalLab/terra/nanoporeTerra/xpore_modRate_DRACH_mettl3kd.bed", sep="\t",col.names = F)


xptable.gr %>% as.data.frame() %>% filter(mod_rate_CTRLSH.REP1 > 0 & grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T)) %>% dplyr::select(kmer,mod_rate_CTRLSH.REP1) %>% distinct() %>% dplyr::group_by(kmer) %>% ggplot(aes(fct_reorder(kmer,mod_rate_CTRLSH.REP1),mod_rate_CTRLSH.REP1)) + geom_boxplot(alpha=0.5,outlier.alpha = 0) + geom_point(alpha=0.2,size=0.5, position = position_dodge2(width = 0.3)) + theme_bw() + labs(x="k-mer (DRACH)", y="Modification rate Ctrl-sh")


modRatesAtChrEnds_csh %>% as.data.frame() %>% filter(mod_rate_CTRLSH.REP1 > 0 & grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T)) %>% dplyr::select(kmer,mod_rate_CTRLSH.REP1) %>% distinct() %>% dplyr::group_by(kmer) %>% ggplot(aes(fct_reorder(kmer,mod_rate_CTRLSH.REP1),mod_rate_CTRLSH.REP1)) + geom_boxplot(alpha=0.5,outlier.alpha = 0.7) + geom_point(alpha=0.7,size=1) + theme_bw() + labs(x="k-mer (DRACH)", y="Modification rate Ctrl-sh")

xptable.gr %>% as.data.frame() %>% filter(mod_rate_METTL3KD.REP1 > 0 & grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T)) %>% dplyr::select(kmer,mod_rate_METTL3KD.REP1) %>% distinct() %>% dplyr::group_by(kmer) %>% ggplot(aes(fct_reorder(kmer,mod_rate_METTL3KD.REP1),mod_rate_METTL3KD.REP1)) + geom_boxplot(alpha=0.5,outlier.alpha = 0) + geom_point(alpha=0.2,size=0.5, position = position_dodge2(width = 0.3)) + theme_bw() + labs(x="k-mer (DRACH)", y="Modification rate METTL3-KD")

#NNANN
modRatesAtChrEnds_mettl3kd %>% as.data.frame() %>% filter(mod_rate_METTL3KD.REP1 > 0 & grepl("..AC.",kmer,perl = T)) %>% dplyr::select(kmer,mod_rate_METTL3KD.REP1) %>% distinct() %>% dplyr::group_by(kmer) %>% ggplot(aes(fct_reorder(kmer,mod_rate_METTL3KD.REP1),mod_rate_METTL3KD.REP1)) + geom_boxplot(alpha=0.5,outlier.alpha = 0.7) + geom_point(alpha=0.7,size=1) + theme_bw() + labs(x="k-mer (NNACN)", y="Modification rate METTL3-KD") + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))


modRatesAtChrEnds_csh %>% as.data.frame() %>% filter(mod_rate_CTRLSH.REP1 > 0 & grepl("..AC.",kmer,perl = T)) %>% dplyr::select(kmer,mod_rate_CTRLSH.REP1) %>% distinct() %>% dplyr::group_by(kmer) %>% ggplot(aes(fct_reorder(kmer,mod_rate_CTRLSH.REP1),mod_rate_CTRLSH.REP1)) + geom_boxplot(alpha=0.5,outlier.alpha = 0.7) + geom_point(alpha=0.7,size=1) + theme_bw() + labs(x="k-mer (NNACN)", y="Modification rate Ctrl-sh") + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))

xptable.gr %>% as.data.frame() %>% filter(mod_rate_CTRLSH.REP1 > 0 & grepl("..AC.",kmer,perl = T)) %>% dplyr::select(kmer,mod_rate_CTRLSH.REP1) %>% distinct() %>% dplyr::group_by(kmer) %>% ggplot(aes(fct_reorder(kmer,mod_rate_CTRLSH.REP1),mod_rate_CTRLSH.REP1)) + geom_boxplot(alpha=0.5,outlier.alpha = 0) + geom_point(alpha=0.2,size=0.5, position = position_dodge2(width = 0.3)) + theme_bw() + labs(x="k-mer (NNACN)", y="Modification rate Ctrl-sh") + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))


xptable.gr %>% as.data.frame() %>% filter(mod_rate_METTL3KD.REP1 > 0 & grepl("..AC.",kmer,perl = T)) %>% dplyr::select(kmer,mod_rate_METTL3KD.REP1) %>% distinct() %>% dplyr::group_by(kmer) %>% ggplot(aes(fct_reorder(kmer,mod_rate_METTL3KD.REP1),mod_rate_METTL3KD.REP1)) + geom_boxplot(alpha=0.5,outlier.alpha = 0) + geom_point(alpha=0.2,size=0.5, position = position_dodge2(width = 0.3)) + theme_bw() + labs(x="k-mer (NNACN)", y="Modification rate METTL3-KD") + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))

#| G(A|G)A(T|G)(T|G) non-canonical motif Xu et al 2021:

modRatesAtChrEnds_csh %>% as.data.frame() %>% filter(mod_rate_CTRLSH.REP1 > 0.01 & grepl("((TAA(A|G)(A|G))| (G(A|G)A(T|G)(T|G)))",kmer,perl = T)) %>% dplyr::select(kmer,mod_rate_CTRLSH.REP1) %>% distinct() %>% dplyr::group_by(kmer) %>% ggplot(aes(fct_reorder(kmer,mod_rate_CTRLSH.REP1),mod_rate_CTRLSH.REP1)) + geom_boxplot(alpha=0.5,outlier.alpha = 0.7) + geom_point(alpha=0.7,size=1) + theme_bw() + labs(x="k-mer (TAA(A|G)(A|G) | G(A|G)A(T|G)(T|G) non-canonical motif)", y="Modification rate CTRLSH") + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))


modRatesAtChrEnds_mettl3kd %>% as.data.frame() %>% filter(mod_rate_METTL3KD.REP1 > 0.01 & grepl("((TAA(A|G)(A|G))| (G(A|G)A(T|G)(T|G)))",kmer,perl = T)) %>% dplyr::select(kmer,mod_rate_METTL3KD.REP1) %>% distinct() %>% dplyr::group_by(kmer) %>% ggplot(aes(fct_reorder(kmer,mod_rate_METTL3KD.REP1),mod_rate_METTL3KD.REP1)) + geom_boxplot(alpha=0.5,outlier.alpha = 0.7) + geom_point(alpha=0.7,size=1) + theme_bw() + labs(x="k-mer (TAA(A|G)(A|G) | G(A|G)A(T|G)(T|G) non-canonical motif)", y="Modification rate METTL3-KD") + theme(axis.text.x = element_text(angle = 45,vjust = 1,hjust = 1))

# Differential modification rates at chr ends, NNACN motif:

modRatesAtChrEnds_csh %>% as.data.frame() %>% arrange(desc(mod_rate_CTRLSH.REP1)) %>% filter(pval_CTRLSH_vs_METTL3KD < 0.05 & grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T)) %>% pivot_longer(names_to = "sample", values_to = "mod_rate",c("mod_rate_CTRLSH.REP1","mod_rate_METTL3KD.REP1")) %>% mutate(sample=gsub("(mod_rate_|.REP1)","",sample,perl = T)) %>% ggplot(aes(kmer,mod_rate,color=sample)) + geom_boxplot(alpha=0.5) + geom_point(aes(color=sample),size=0.75,alpha=0.75, position = position_dodge(0.8)) + theme_bw() + labs(x="k-mer",y="Modification rate", title = "xPore Modification rates Ctrl-sh and METTL3-KD \n chr ends 50kb, NNACN")


# Genome-wide Differential modification rates at chr ends, all motifs:

topSigkmer_drach <- xptable.gr %>% filter(grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T)) %>% as.data.frame() %>% arrange(desc(-log10(pval_CTRLSH_vs_METTL3KD))) %>% head(100) %>% dplyr::select(kmer)

xptable.gr %>% as.data.frame() %>%
  filter(grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T)) %>%
  mutate(label=if_else(-log10(pval_CTRLSH_vs_METTL3KD) > 25 & abs(diff_mod_rate_CTRLSH_vs_METTL3KD) > 0.7,kmer,"")) %>% 
  ggplot(aes(diff_mod_rate_CTRLSH_vs_METTL3KD,-log10(pval_CTRLSH_vs_METTL3KD),label=label)) + 
  geom_text_repel(size=2) +
  geom_point(size=0.05, alpha=0.3) + theme_bw() + labs(color="p<0.05",x="Differential modification rate", y="-log10(p-value)", title = "Diff. modification rate Ctrl-sh vs METTL3-KD \n DRACH motif") + ylim(0,100)

#%>% filter(mod_rate_CTRLSH.REP1 > 0.5 & coverage_CTRLSH.REP1 > 5) %>% dplyr::select(kmer,mod_rate_CTRLSH.REP1,position) %>% group_by(kmer) %>% arrange(desc(mod_rate_CTRLSH.REP1)) %>% top_n(1) 
#%>% mutate(color=if_else(grepl("(A|G|T)(A|G)AC(A|C|T)",kmer),"red",if_else(grepl("(GGATT|GAAGG|TAAAG|TAAGA)"))))

 #%>% ggplot(aes(reorder(kmer,dplyr::desc(mod_rate_CTRLSH.REP1)),mod_rate_CTRLSH.REP1)) + geom_point() + theme(axis.text.x = element_text(angle=45,vjust = 1,hjust = 1,size = 7))


peaksCsh <- rtracklayer::import("~/MondalLab/terra/gathered_peaks/Csh_noTreat_ext75_peaks_all.narrowPeak")

tmp <- join_overlap_inner(xptable.gr ,peaksCsh)



xptable_genwide_t2t <- fread("~/MondalLab/terra/nanoporeTerra/xpore_diffmod/diffmod.table", sep = ",",header = T)

tmp <- xptable_genwide_t2t %>% filter( diff_mod_rate_CTRLSH_vs_METTL3KD >  0, 
                                            pval_CTRLSH_vs_METTL3KD < 0.05,
                                            #`mod_rate_CTRLSH-REP1` > 0.5,
                                            grepl("..A..",kmer)) 


left_join(tmp %>% dplyr::select(id,position,kmer, `mod_rate_CTRLSH-REP1`), genePred, by=c("id"="gene_id")) %>%
  filter(!is.na(gene_start)) %>% 
  mutate(start=position-10,end=position+10,range=paste0(chr,":",start,"-",end)) %>% 
  dplyr::select(chr,start,end,kmer,`mod_rate_CTRLSH-REP1`, gene_strand) %>% 
  fwrite(.,"~/MondalLab/terra/nanoporeTerra/predicted_sites_xpore_lists/xpore_t2t_genomeWide_diffMod_pval0.05_NNANN_10bpAround.bed",sep = "\t",col.names = F)


tmp1 <- left_join(tmp %>% dplyr::select(id,position,kmer, `mod_rate_CTRLSH-REP1`), genePred, by=c("id"="gene_id")) %>%
  filter(!is.na(gene_start)) %>% 
  mutate(start=position-10,end=position+10,range=paste0(chr,":",start,"-",end)) %>% 
  dplyr::select(chr,start,end,kmer,`mod_rate_CTRLSH-REP1`, gene_strand) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)

join_overlap_inner(tmp1,chrEndCoords_30kb)


top100 <- tmp %>% mutate(row=row_number()) %>% arrange(desc(-log10(pval_CTRLSH_vs_METTL3KD))) %>% head(100) %>% dplyr::select(row) %>% unlist()

top100_kmers <- tmp %>% mutate(row=row_number()) %>% arrange(desc(-log10(pval_CTRLSH_vs_METTL3KD))) %>% head(100) %>% dplyr::select(kmer) %>% unlist()

all_kmers <- tmp %>% dplyr::select(kmer) %>% unlist()

library(ggseqlogo)

ggseqlogo(all_kmers, method="prob")

tmp %>% mutate(row=row_number()) %>% mutate(label=if_else(row %in% top100,kmer,NULL)) %>% 
  ggplot(aes(diff_mod_rate_CTRLSH_vs_METTL3KD,-log10(pval_CTRLSH_vs_METTL3KD), label=label)) +
  geom_point(size=0.01, alpha=0.5) + 
  geom_text_repel(size=2, max.overlaps = Inf,  segment.alpha=0.3, segment.size=0.3) +
  scale_y_continuous(limits =c(0,NA)) + theme_prism(base_size = 10) +
  labs(title="Differentially modified sites, xPore T2T", x="Differential modification rate \n Control-sh vs METTL3-KD", y="-log10(p-value)") -> p

ggsave("~/MondalLab/terra/plots/figure_compilation/differentially_modified_sites_xPore_T2T_genomewide_diffmodRateAbove0_pval0.05_NNANN.pdf", plot = p)

#((A|G|T)(A|G)AC(A|C|T))| ((A|G)AC(A|C|T)(A|G|T)) | (AC(A|C|T)(A|G|T)(A|G)) | (C(A|C|T)(A|G|T)(A|G)A) | ((A|C|T)(A|G|T)(A|G)AC)

tmp %>%
  #filter(diff_mod_rate_CTRLSH_vs_METTL3KD > 0.9) %>% 
  arrange(desc(-log10(pval_CTRLSH_vs_METTL3KD))) %>% 
  #head(100) %>%
  mutate(drach=if_else(grepl("..AC.",kmer,perl = T),"DRACH","Non-DRACH")) %>% #(A|G|T)(A|G)AC(A|C|T)
 group_by(drach) %>% 
  count() %>%
  ggplot(aes(fct_reorder(drach,rev(c("DRACH","Non-DRACH"))),n, fill=drach)) + geom_col(alpha=0.35) +theme_prism(axis_text_angle = 0, base_size = 10) +
  labs(y="Number of differentially modified sites", x="k-mer pattern") + theme(legend.position = "none") + scale_fill_manual(values = c("steelblue","black")) + coord_flip()
  

tmp %>% 
  #filter(diff_mod_rate_CTRLSH_vs_METTL3KD > 0.9) %>% 
  arrange(desc(-log10(pval_CTRLSH_vs_METTL3KD))) %>% 
  #head(100) %>%
  mutate(drach=if_else(grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T),TRUE,FALSE)) %>% filter(drach==TRUE) %>% dplyr::select(kmer) %>% ggseqlogo(.,method="prob") + labs(title="DRACH")
  
  # 
  # mutate(medianModRate=mean(diff_mod_rate_CTRLSH_vs_METTL3KD)) %>% 
  # dplyr::select(kmer,medianModRate) %>% ungroup() %>% distinct() %>%
  # arrange(desc(medianModRate)) %>% 
  # mutate(drach=if_else(grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T),TRUE,FALSE)) %>% 
  # ggplot(aes(reorder(kmer, desc(medianModRate)), medianModRate, fill=drach)) + geom_col() + theme_prism(axis_text_angle = 45, base_size = 10)

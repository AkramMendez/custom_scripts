library(tidyverse)
library(data.table)
library(plyranges)
library(rtracklayer)


chrEndCoords <- rtracklayer::import("~/gitrepos/bioinfopipes/lab/CHM13_hg38chrY_chromEndCoordinates_50kb.genome",format="BED")

chrEndCoords <- chrEndCoords %>% as.data.frame() %>% mutate(strand=if_else(end > 50000,"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)

# Load xPore diffmod table resutls, for a description of each column check: https://xpore.readthedocs.io/en/latest/outputtable.html#outputtable
xptable <- fread("~/MondalLab/terra/nanoporeTerra/xpore_diffmod/xpore_diffmod_t2tTxm_chrEnds50kb/diffmod.table", sep = ",",header = T)

#xptable %>% filter(`mod_rate_CTRLSH-REP1` > 0.8, conf_mu_mod > 0.2, mod_assignment=="higher",diff_mod_rate_CTRLSH_vs_METTL3KD > 0.5 )

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



#(G|A|T)GT(T|C)(T|C)
modRatesAtChrEnds_csh <- join_overlap_inner_directed(xptable.gr,chrEndCoords)

modRatesAtChrEnds_csh

#modRatesAtChrEnds_csh 
xptable.gr %>% as.data.frame() %>% filter(diff_mod_rate_CTRLSH_vs_METTL3KD > 0,mod_rate_CTRLSH.REP1 > 0.9, mod_assignment=="higher", coverage_CTRLSH.REP1 > 10) %>% arrange(desc(z_score_CTRLSH_vs_METTL3KD)) %>% ggplot(aes(diff_mod_rate_CTRLSH_vs_METTL3KD,mod_rate_CTRLSH.REP1, label=kmer)) + geom_point(size=1) + geom_text_repel(size=3, max.overlaps = Inf) + scale_x_continuous(expand = expansion(mult = c(0,1.5))) +scale_y_continuous(expand = expansion(mult = c(0,1.5)))

# Boxplot, modification rates
xptable.gr %>% as.data.frame() %>% filter( diff_mod_rate_CTRLSH_vs_METTL3KD > 0, grepl("..A..",kmer,perl = T), pval_CTRLSH_vs_METTL3KD < 0.05 ) %>% pivot_longer(names_to = "sample",values_to = "mod_rate", matches("\\bmod_rate_")) %>% dplyr::select(sample,mod_rate, kmer,seqnames,start,end,strand) %>% ggplot(aes(sample,mod_rate)) + geom_quasirandom(width = 0.25) + geom_boxplot(alpha=0, width = 0.5) + theme_prism()

#xptable.gr %>% as.data.frame() %>% filter(mod_rate_CTRLSH.REP1 > 0 & grepl("..A..",kmer,perl = T)) %>% dplyr::select(seqnames,start,end,kmer,mod_rate_CTRLSH.REP1,strand) %>% fwrite(.,"~/MondalLab/terra/nanoporeTerra/xpore_modRateAbove0.9_ctrlsh_withSecondary_chrEnds50kb.bed", sep="\t",col.names = F)
#mod_assignment=="higher"
xptable.gr %>% as.data.frame() %>% filter( diff_mod_rate_CTRLSH_vs_METTL3KD > 0.25 & grepl("..A..",kmer,perl = T) ) %>% dplyr::select(seqnames,start,end,kmer,mod_rate_CTRLSH.REP1,strand) %>% fwrite(.,"~/MondalLab/terra/nanoporeTerra/xpore_modRateAbove0.9_ctrlsh_withSecondary_chrEnds50kb.bed", sep="\t",col.names = F)

xptable.gr %>% as.data.frame() %>% filter( diff_mod_rate_CTRLSH_vs_METTL3KD > 0) %>% dplyr::select(seqnames,start,end,kmer,mod_rate_METTL3KD.REP1,strand) %>% fwrite(.,"~/MondalLab/terra/nanoporeTerra/xpore_modRateAbove0.9_mettl3kd_withSecondary_chrEnds50kb.bed", sep="\t",col.names = F)


xptable.gr %>% as.data.frame() %>% filter(mod_rate_CTRLSH.REP1 > 0.9, coverage_CTRLSH.REP1 > 5 & grepl("..A..",kmer,perl = T) & diff_mod_rate_CTRLSH_vs_METTL3KD > 0) %>% dplyr::select(seqnames,start,end,kmer,mod_rate_CTRLSH.REP1,strand) %>% fwrite(.,"~/MondalLab/terra/nanoporeTerra/xpore_modRateAbove0.9_NNANN_ctrlsh_withSecondary_chrEnds50kb.bed", sep="\t",col.names = F)

xptable.gr %>% as.data.frame() %>% filter(mod_rate_METTL3KD.REP1 > 0.9, coverage_METTL3KD.REP1 > 5 & grepl("..A..",kmer,perl = T) & diff_mod_rate_CTRLSH_vs_METTL3KD > 0) %>% dplyr::select(seqnames,start,end,kmer,mod_rate_METTL3KD.REP1,strand) %>% fwrite(.,"~/MondalLab/terra/nanoporeTerra/xpore_modRateAbove0.9_NNANN_mettl3kd_withSecondary_chrEnds50kb.bed", sep="\t",col.names = F)


xptable.gr %>% as.data.frame() %>% filter(mod_rate_CTRLSH.REP1 > 0.5, coverage_CTRLSH.REP1 > 5 & grepl("..A..",kmer,perl = T) & diff_mod_rate_CTRLSH_vs_METTL3KD > 0) %>% dplyr::select(seqnames,start,end,kmer,mod_rate_CTRLSH.REP1,strand) %>% fwrite(.,"~/MondalLab/terra/nanoporeTerra/xpore_modRateAbove0.5_NNANN_ctrlsh_withSecondary_chrEnds50kb.bed", sep="\t",col.names = F)

xptable.gr %>% as.data.frame() %>% filter(mod_rate_METTL3KD.REP1 > 0.5, coverage_METTL3KD.REP1 > 5 & grepl("..A..",kmer,perl = T) & diff_mod_rate_CTRLSH_vs_METTL3KD > 0) %>% dplyr::select(seqnames,start,end,kmer,mod_rate_METTL3KD.REP1,strand) %>% fwrite(.,"~/MondalLab/terra/nanoporeTerra/xpore_modRateAbove0.5_NNANN_mettl3kd_withSecondary_chrEnds50kb.bed", sep="\t",col.names = F)






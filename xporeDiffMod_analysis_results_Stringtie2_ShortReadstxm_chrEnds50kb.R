library(tidyverse)
library(data.table)
library(plyranges)
library(rtracklayer)


chrEndCoords <- rtracklayer::import("~/gitrepos/bioinfopipes/lab/CHM13_hg38chrY_chromEndCoordinates_50kb.genome",format="BED")

chrEndCoords <- chrEndCoords %>% as.data.frame() %>% mutate(strand=if_else(end > 50000,"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)


# Load xPore diffmod table resutls, for a description of each column check: https://xpore.readthedocs.io/en/latest/outputtable.html#outputtable
xptable_str_shorReads <- fread("~/MondalLab/terra/nanoporeTerra/xpore_diffmod/xpore_diffmod_stringtieShortReadsTxm_chrEnds50kb/diffmod.table", sep = ",",header = T) %>% mutate(id=gsub("_",".",id))

#xptable %>% filter(`mod_rate_CTRLSH-REP1` > 0.8, conf_mu_mod > 0.2, mod_assignment=="higher",diff_mod_rate_CTRLSH_vs_METTL3KD > 0.5 )

## Convert Stringtie GTF to genePred for extracting transcript coordinates (on terminal) :
#gtfToGenePred -genePredExt ctrlsh_withSecondary_sp.chm13.genome.14.sorted_txm_Stringtie2.2.1_notUCSCformat_underscore.gtf ctrlsh_withSecondary_sp.chm13.genome.14.sorted_txm_Stringtie2.2.1_notUCSCformat_underscore_gffread.genePred
#awk 'BEGIN{FS=OFS="\t"}{if($7=="."){print $1,$2,$3,$4,$5,$6,"+",$7,$8,$9}else{print}}' csh_merged_input.nodup.sorted_txm_Stringtie2.2.1.gtf > csh_merged_input.nodup.sorted_txm_Stringtie2.2.1_strandInfo.gtf
#gff_str <- rtracklayer::import("/home/akram/gitrepos/bioinfopipes/lab/ont_scripts/csh_merged_input.nodup.sorted_txm_Stringtie2.2.1.gtf")
gff_str_shortReads <- makeTxDbFromGFF("/home/akram/gitrepos/bioinfopipes/lab/ont_scripts/csh_merged_input.nodup.sorted_txm_Stringtie2.2.1_strandInfo.gtf",format = "gtf")
#gff_str <- gff_str %>% filter(seqnames!="U00096.2")

#gff_str <- gff_str %>% as.data.frame() %>% dplyr::select(seqnames,start,end,width,strand,gene_id) %>% distinct()
genePred_str <- fread("/home/akram/gitrepos/bioinfopipes/lab/ont_scripts/csh_merged_input.nodup.sorted_txm_Stringtie2.2.1.genePred",sep = "\t", header = F)[,c(1:5,12)]

colnames(genePred_str) <- c("tx_id","chr","gene_strand","start","end","gene_id")
genePred_str <- genePred_str %>% filter(chr!="U00096.2")

genePred_str <- genePred_str %>% 
  group_by(gene_id) %>% 
  mutate(gene_start=min(start),gene_end=max(end)) %>%
  dplyr::select(chr,gene_start,gene_end,gene_strand,gene_id,tx_id) %>%
  distinct() %>% 
  ungroup()

genePred_str <- genePred_str %>% mutate(gene_width=gene_end-gene_start)

xptable_str_shorReads <- left_join(xptable_str_shorReads,genePred_str,by=c("id"="tx_id"))
#left_join(xptable_str,genePred_str,by=c("id"="gene_id"))
#"(A|G|T)(A|G)AC(A|C|T)"
xptable_str_shorReads <-
  xptable_str_shorReads %>% mutate(drach = if_else(grepl("(A|G|T)(A|G)AC(A|C|T)", kmer, perl = T), TRUE, FALSE))

xptable_str_shorReads <-
  xptable_str_shorReads %>% filter(!is.na(gene_start))

xptable_str_shorReads <- xptable_str_shorReads %>% mutate(
  start = if_else(
    gene_strand == "+",
    (gene_start + position) - 3,
    (gene_end - position) - 3
  ),
  end = if_else(
    gene_strand == "+",
    (gene_start + position) + 3,
    (gene_end - position) + 3
  )
) %>% mutate(width = end - start)
#mutate(start=position-3,end=position+2) %>% mutate(width=end-start)

#%>% dplyr::select(chr,gene_start,gene_end,gene_width,gene_strand,colnames(xptable)[1:19],drach)
xptable_str_shorReads.gr <- xptable_str_shorReads %>% 
  #dplyr::select(chr,start,end,width,gene_strand,colnames(xptable_str)[1:19],drach) %>% 
  dplyr::select(chr,gene_start,gene_end,gene_width,gene_strand,colnames(xptable)[1:19],drach) %>%
  rename(start=gene_start,end=gene_end) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T,seqnames.field = "chr",strand.field = "gene_strand",end.field = "end",start.field = "start")




join_overlap_inner_directed(xptable_str_shorReads.gr,chrEndCoords)


xptable_str_shorReads.gr %>% as.data.frame() %>% dplyr::select(kmer,mod_rate_CTRLSH.REP1,position) %>% filter(grepl("((..AC.)|GGATT|GAAGG|TAAAG|TAAGA)",kmer,perl = T)) %>% arrange(dply::desc(mod_rate_CTRLSH.REP1)) %>% head(100)

xptable_str_shorReads.gr %>% as.data.frame() %>% mutate(arm=if_else(end < 100000,"P","Q"),chr=paste0("chr",seqnames,"_",arm)) %>% dplyr::count(chr)



#(G|A|T)GT(T|C)(T|C)
modRatesAtChrEnds_str_shorReads_csh <- join_overlap_inner_directed(xptable_str_shorReads.gr,chrEndCoords)

modRatesAtChrEnds_str_shorReads_csh

#modRatesAtChrEnds_csh (A|G|T)(A|G)AC(A|C|T)
xptable_str_shorReads.gr %>% as.data.frame() %>%
  filter( diff_mod_rate_CTRLSH_vs_METTL3KD > 0, grepl("..A..",kmer,perl = T) , pval_CTRLSH_vs_METTL3KD < 0.05 , mod_rate_CTRLSH.REP1 > 0.25) %>% #t.test < 0.05
  arrange(desc(z_score_CTRLSH_vs_METTL3KD)) %>% ggplot(aes(diff_mod_rate_CTRLSH_vs_METTL3KD,mod_rate_CTRLSH.REP1, label=kmer)) + geom_point(size=1) + geom_text_repel(size=3, max.overlaps = Inf) + scale_x_continuous(expand = expansion(mult = c(0,.1))) +scale_y_continuous(expand = expansion(mult = c(0,.1)))

# Boxplot, modification rates
xptable_str_shorReads.gr %>% as.data.frame() %>%
  filter( diff_mod_rate_CTRLSH_vs_METTL3KD > 0, grepl("..A..",kmer,perl = T), pval_CTRLSH_vs_METTL3KD < 0.05 , mod_rate_CTRLSH.REP1 > 0.25, t.test < 0.05) %>%
  mutate(arm=if_else(end < 100000,"P","Q"),chr=paste0("chr",seqnames,"_",arm)) %>% 
  dplyr::select(chr,id,kmer,mod_rate_CTRLSH.REP1,mod_rate_METTL3KD.REP1,diff_mod_rate_CTRLSH_vs_METTL3KD) %>% count(chr)

df <- xptable_str_shorReads.gr %>% as.data.frame() %>% 
  filter( diff_mod_rate_CTRLSH_vs_METTL3KD > 0, grepl("..A..",kmer,perl = T), pval_CTRLSH_vs_METTL3KD < 0.05 , mod_rate_CTRLSH.REP1 > 0.25, t.test < 0.05) %>% 
  pivot_longer(names_to = "sample",values_to = "mod_rate", matches("\\bmod_rate_")) %>% 
  mutate(sample=gsub("mod_rate_","",sample)) %>% 
  dplyr::select(sample,mod_rate, kmer,seqnames,start,end,strand) 

df$sample <- factor(df$sample,c("CTRLSH.REP1","METTL3KD.REP1"))

df %>% ggpaired(.,x="sample",y="mod_rate", line.color=NA, line.size = 0.4,point.size =1) + theme_prism(axis_text_angle = 45) + scale_y_continuous(expand = expansion(mult = c(0.01,0.1))) + stat_compare_means(ref.group = "CTRLSH.REP1",paired = T,label = "p.signif") + labs(x="",y="Modification rate", title="xPore Stringtie2 transcripts")

# Compare to results from xPore on annotated T2T transcripts
xptable.gr %>% as.data.frame() %>%
  filter( diff_mod_rate_CTRLSH_vs_METTL3KD > 0, grepl("..A..",kmer,perl = T), pval_CTRLSH_vs_METTL3KD < 0.05 , mod_rate_CTRLSH.REP1 > 0.25 ) %>% #t.test < 0.05
  mutate(arm=if_else(end < 100000,"P","Q"),chr=paste0(seqnames,"_",arm)) %>% 
  dplyr::select(chr,id,kmer,mod_rate_CTRLSH.REP1,mod_rate_METTL3KD.REP1,diff_mod_rate_CTRLSH_vs_METTL3KD) %>% count(chr)

df_xpore_t2t <- xptable.gr %>% as.data.frame() %>% 
  filter( diff_mod_rate_CTRLSH_vs_METTL3KD > 0, grepl("..A..",kmer,perl = T), pval_CTRLSH_vs_METTL3KD < 0.05 , mod_rate_CTRLSH.REP1 > 0.25) %>% 
  pivot_longer(names_to = "sample",values_to = "mod_rate", matches("\\bmod_rate_")) %>% 
  mutate(sample=gsub("mod_rate_","",sample), seqnames=as.factor(gsub("chr","",seqnames))) %>% 
  dplyr::select(sample,mod_rate, kmer,seqnames,start,end,strand) 

df_xpore_t2t$sample <- factor(df_xpore_t2t$sample,c("CTRLSH.REP1","METTL3KD.REP1"))

df_xpore_t2t %>% 
  ggpaired(.,x="sample",y="mod_rate", line.color=NA, line.size = 0.4,point.size =1) + theme_prism(axis_text_angle = 45) + scale_y_continuous(expand = expansion(mult = c(0.01,0.1))) + 
  stat_compare_means(ref.group = "CTRLSH.REP1",paired = T,label = "p.signif") + labs(x="",y="Modification rate", title="xPore T2T transcripts")

# Select predicted sites from xPore results using Stringtie short-reads assembly
modified_stringtie_shortReads <- xptable_str_shorReads.gr %>% as.data.frame() %>% 
  filter( diff_mod_rate_CTRLSH_vs_METTL3KD > 0, grepl("..A..",kmer,perl = T), pval_CTRLSH_vs_METTL3KD < 0.05 , mod_rate_CTRLSH.REP1 > 0.25, t.test < 0.05) %>%
  mutate(arm=if_else(end < 100000,"P","Q"),chr=paste0("chr",seqnames,"_",arm)) %>% dplyr::select(chr) %>% distinct() %>% unlist()

# Select predicted sites from xPore results using annotated T2T transcripts
modified_t2t <- xptable.gr %>% as.data.frame() %>% 
  filter( diff_mod_rate_CTRLSH_vs_METTL3KD > 0, grepl("..A..",kmer,perl = T), pval_CTRLSH_vs_METTL3KD < 0.05 , mod_rate_CTRLSH.REP1 > 0.25, t.test < 0.05) %>%
  mutate(arm=if_else(end < 100000,"P","Q"),chr=paste0(seqnames,"_",arm)) %>% dplyr::select(chr) %>% distinct() %>% unlist()

firing_chroms <- paste0("chr",c("11_P","14_P","9_P","15_Q","3_Q","20_Q","12_P","16_P","19_P","21_Q"))

join_overlap_inner(peaksCsh,chrEndCoords) %>% as.data.frame() %>% mutate(arm=if_else(end > 100000,"Q","P")) %>% mutate(label=paste0(seqnames,"_",arm)) %>% dplyr::select(label) %>% distinct() %>% unlist() -> arms_chrEnd50kb_m6aRIP_ctrlsh_u2os


firing_modified_stringtie_shorReads_t2t_m6arip <- intersect(intersect(modified_stringtie_shorReads,modified_t2t),arms_chrEnd50kb_m6aRIP_ctrlsh_u2os)#,arms_chrEnd50kb_m6aRIP_ctrlsh_u2os),firing_chroms)



# Convert trasncriptomic to genomic coordinates
colnames(xptable_str_shorReads) <- gsub("-",".",colnames(xptable_str_shorReads))

coords_shortReads <- xptable_str_shorReads %>% #Filter results to select significant modifications, NNANN context
  filter( diff_mod_rate_CTRLSH_vs_METTL3KD > 0, grepl("..A..",kmer,perl = T), pval_CTRLSH_vs_METTL3KD < 0.05) %>%
  dplyr::select(id,position,kmer,diff_mod_rate_CTRLSH_vs_METTL3KD,pval_CTRLSH_vs_METTL3KD,mod_rate_CTRLSH.REP1,mod_rate_METTL3KD.REP1) %>% 
  mutate(seqnames=id,start=position-3,end=position+3,width=end-start, position_in_tx=position) %>%
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

transcripts_str_shorReads <- transcripts(gff_str_shortReads)
exons_str_shorReads <- GenomicFeatures::exons(gff_str_shortReads,columns=c("exon_id","tx_name","gene_id"))
exons_str_shorReads <- exons_str_shorReads %>% as.data.frame() %>% separate_rows(.,tx_name,gene_id,sep = " ,")
names(transcripts_str_shorReads) <- transcripts_str_shorReads$tx_name


mapping_shorReads <- mapFromTranscripts(coords_shortReads, transcripts_str_shorReads,ignore.strand=FALSE) %>% as.data.frame() %>% mutate(tx_id=coords_shortReads$id, kmer=coords_shortReads$kmer, position_in_tx=coords_shortReads$position_in_tx) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)


#Merge exon information with merged trasncript/xPore predictions mapping:

mapping_and_exons_shortReads <- left_join(
  mapping_shorReads %>% as.data.frame(),
  exons_str_shorReads %>% as.data.frame() %>% rename_with( ~ paste0(.x, "_exon")),
  by = c("tx_id" = "tx_name_exon")
)



translateToGenome <- function(mapping_and_exons){
  
  neg_strand <- mapping_and_exons %>% filter(strand=="-")
  pos_strand <- mapping_and_exons %>% filter(strand=="+")
  
  #Negative strand
  df_neg <- neg_strand %>%
    arrange(kmer) %>% 
    group_by(kmer) %>% 
    group_by(tx_id,position_in_tx,kmer) %>%
    arrange(kmer,desc(exon_id_exon)) %>% 
    mutate(tx_len=sum(width_exon), 
           until=cumsum(width_exon), 
           which_exon=if_else(position_in_tx <= until,TRUE,FALSE),
           last_until=if_else(position_in_tx < until & all(which_exon==TRUE),until[nth(which(which_exon==TRUE),1L)],until[nth(which(which_exon==FALSE),-1L)]),
           last_end=end_exon[nth(which(which_exon==TRUE),1L)],
           dist_from_start=if_else(position_in_tx < until & all(which_exon==TRUE),position_in_tx,position_in_tx-nth(last_until,-2L)),
           coord=last_end-dist_from_start,
           chr=paste0("chr",seqnames),
           coord_start=coord-3,
           coord_end=coord+2,
           range=paste0("chr",seqnames,":",coord-3,"-",coord+2)) %>%
    ungroup() %>% 
    as.data.frame()
  
  #Positive strand
  df_pos <- pos_strand %>%
    arrange(kmer) %>% 
    group_by(kmer) %>% 
    group_by(tx_id,position_in_tx) %>%
    #arrange(kmer,desc(exon_id_exon)) %>% 
    mutate(tx_len=sum(width_exon), 
           until=cumsum(width_exon), 
           which_exon=if_else(position_in_tx <= until,TRUE,FALSE),
           last_until=if_else(position_in_tx < until & all(which_exon==TRUE),until[nth(which(which_exon==TRUE),1L)],until[nth(which(which_exon==FALSE),-1L)]),
           last_end=start_exon[nth(which(which_exon==TRUE),1L)],
           dist_from_start=if_else(position_in_tx < until & all(which_exon==TRUE),position_in_tx,position_in_tx-nth(last_until,-2L)),
           coord=last_end+dist_from_start,
           chr=paste0("chr",seqnames),
           coord_start=coord-3,
           coord_end=coord+2,
           range=paste0("chr",seqnames,":",coord-3,"-",coord+2)) %>%
    ungroup() %>% 
    as.data.frame()  
  
  df <- rbind(df_pos,df_neg)
  
  df <- df %>% dplyr::select(chr,coord_start,coord_end,kmer,position_in_tx,strand,tx_id,range) %>% distinct()
  
  return(df)
  
}


passing_sites_genomicCoors_str_shorReads <- translateToGenome(mapping_and_exons_shortReads %>% mutate(seqnames=gsub("chr","",seqnames))) %>% mutate(strand=".") %>% dplyr::select(chr,coord_start,coord_end,kmer,position_in_tx,strand) 

passing_sites_genomicCoors_str_shorReads %>% fwrite("~/MondalLab/terra/nanoporeTerra/xpore_diffmod/xpore_diffmod_stringtieShortReadsTxm_chrEnds50kb/coordinates_tx_to_genome_stringtie_shortReads_xpore_chrEnds50kb.bed", col.names = F,sep="\t")

# Save predicted sites to excel table


# Merge with expression CPM values (object 'df' from the "counts_coverage_cheEnds50kb_CPM.R script)
expression_cpm <- df %>% mutate(group=if_else(which_arm %in% chrmsHigh,"high","low")) %>% dplyr::select(name.y,cpm,which_arm,group) %>% distinct() 

expression_cpm <- rbind(
  expression_cpm,
  df_ONT %>% dplyr::select(name.y,cpm,which_arm,group) %>% distinct()
) %>% pivot_wider(names_from = "name.y", values_from = "cpm") %>% as.data.frame()

colnames(expression_cpm) <- c("chr_arm","expression_class_CtrlshU2OS","cpm_Controlsh_U2OS", "cpm_METTL3_KD_U2OS","cpm_NB_tumor_1","cpm_NB_tumor_2","cpm_Controlsh_SKNF1","cpm_METTL3KD_SKNF1", "cpm_Controlsh_ONT","cpm_METTL3KD_ONT")


left_join(
  left_join(
    translateToGenome(mapping_and_exons_shortReads %>% mutate(seqnames=gsub("chr","",seqnames))),
    coords_shortReads %>% as.data.frame(),
    by = c(
      "tx_id" = "id",
      "position_in_tx" = "position",
      "kmer" = "kmer"
    )
  ) %>% 
    dplyr::select(!matches("(.y$|strand|width|^start$|^end$|seqnames)"))  %>% 
    mutate(label=if_else(coord_end > 100000,paste0(chr,"Q"),paste0(chr,"P"))),
  expression_cpm, by=c("label"="chr_arm")
  
) %>% dplyr::select(!matches("(.y$|strand|width|^start$|^end$|seqnames)")) %>% 
  arrange(pval_CTRLSH_vs_METTL3KD,chr) %>% openxlsx::write.xlsx(.,"~/MondalLab/terra/plots/figure_compilation/tables/predictedSites_xPore_Stringtie_ShortReads_pval0.05_diffModAbove0_NNANN.xlsx", overwrite = T)


# Merge CPM table with Ctrl-sh U2OS chrEnd 50kb peaks from RIP data
join_overlap_inner(
join_overlap_inner(peaksCsh, chrEndCoords),
passing_sites_genomicCoors_str %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)
)


# How many sites are in each chromosome arm? (NNANN)
passing_sites_genomicCoors_str_shorReads %>% mutate(arm=if_else(coord_end > 100000, "Q","P")) %>% mutate(label=paste0(chr,"_",arm)) %>% count(label)

# How many sites are in each chromosome arm? (NNACN) context
passing_sites_genomicCoors_str_shorReads %>% filter(grepl("..AC.",kmer)) %>% mutate(arm=if_else(coord_end > 100000, "Q","P")) %>% mutate(label=paste0(chr,"_",arm)) %>% count(label)
#Intersect with m6A-RIP Control-sh U2OS peaks
join_overlap_inner(passing_sites_genomicCoors_str_shorReads %>% as.data.frame() %>% mutate(start=coord_start-75,end=coord_end+75) %>% dplyr::select(-starts_with("coord")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),
                   peaksCsh,minoverlap = 1)





### Plot modification rates vs differential modification rates NNANN context:

xptable_str_shorReads %>% 
  mutate( label=if_else(diff_mod_rate_CTRLSH_vs_METTL3KD > 0 & grepl("..A..",kmer,perl = T) & pval_CTRLSH_vs_METTL3KD < 0.05 & mod_rate_CTRLSH.REP1 > 0.25,kmer,NULL)) %>% 
  ggplot(aes(diff_mod_rate_CTRLSH_vs_METTL3KD,-log10(pval_CTRLSH_vs_METTL3KD), label=label )) + 
  geom_point(alpha=0.5,size=0.3) +
  geom_text_repel(max.overlaps = Inf, size=3)












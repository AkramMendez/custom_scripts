library(tximport)
library(tximeta)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Guitar)
library(rtracklayer)
library(plyranges)
library(GenomicFeatures)
library(AnnotationDbi)
library(ggprism)
library(ggpubr)

# Custom Functions

readPeaks_xls <- function(xlsPath){
  colnamesPeaks <- c("chr","start","end","length","peak.point.source","pileup","-log10pvalue","fold_enrichment","-log10qvalue","name")
  peaks_gr <- fread(xlsPath, sep="\t", col.names = colnamesPeaks) %>%
    mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
    makeGRangesFromDataFrame(.,keep.extra.columns = T)
  
  return(peaks_gr)
}


myMetagenePlot <- function(df, binwidth){
  df %>% ggplot(aes(x=scaled_summit, color=cond)) + 
    geom_line(stat = "bin", binwidth=binwidth) +
    #geom_density(aes(y=..count..),bw=binwidth) +
    theme_prism(base_size = 10,border = T) +
    scale_y_continuous(expand = expansion(mult=0.01), limits = c(0,NA))-> p
  return(p)
}

overlapPeakAndGenomeCoords<- function(df,genesGTF){
  scaled_df <- join_overlap_inner_directed(df,genesGTF) %>% 
    as.data.frame() %>% 
    mutate(
      summit=start+peak.point.source,
      scaled_start=if_else(
        strand=="+",
        ((start-mg_start)/(mg_end-mg_start))*tot.bins,
        (1-( (start-mg_start)/(mg_end-mg_start) ))*tot.bins),
      scaled_end=if_else(
        strand=="+", 
        (1-(mg_end-end)/(mg_end-mg_start))*tot.bins,
        (((mg_end-end)/(mg_end-mg_start)))*tot.bins),
      scaled_summit=if_else(strand=="+",
                            ((summit-mg_start)/(mg_end-mg_start))*tot.bins,
                            (1-( (summit-mg_start)/(mg_end-mg_start) ))*tot.bins)
    )
  
  return(scaled_df) 
}


overlapPeakAndGenomeCoords2<- function(df,genesGTF, tot.bins){
  
  scaled_df <- join_overlap_inner_within_directed(df,genesGTF,maxgap = 0,minoverlap = 1) %>%
    as.data.frame() %>% 
    mutate(
      summit=start+peak.point.source,
      scaled_start=if_else(
        strand=="+",
        ((start-mg_start)/(mg_end-mg_start))*tot.bins,
        (1-( (start-mg_start)/(mg_end-mg_start) ))*tot.bins),
      scaled_end=if_else(
        strand=="+", 
        (1-(mg_end-end)/(mg_end-mg_start))*tot.bins,
        (((mg_end-end)/(mg_end-mg_start)))*tot.bins),
      scaled_summit=if_else(strand=="+",
                            ((summit-mg_start)/(mg_end-mg_start))*tot.bins,
                            (1-( (summit-mg_start)/(mg_end-mg_start) ))*tot.bins)
    ) 
  
  scaled_df <- scaled_df %>%
    mutate(
      scaled_summit_bin= if_else(
        type=="upstream",
        rescale(scaled_summit, to=c(1,2),from = c(0,tot.bins)),
        if_else(type=="fiveUTR",
                rescale(scaled_summit, to=c(2,3),from = c(0,tot.bins)),
                if_else(type=="cds", 
                        rescale(scaled_summit, to=c(3,6),from = c(0,tot.bins)),
                        if_else(type=="threeUTR",
                                rescale(scaled_summit, to=c(6,7),from = c(0,tot.bins)),
                                if_else(type=="downstream",
                                        rescale(scaled_summit, to=c(7,8),from = c(0,tot.bins)),
                                        0
                                )))))
    )
  
  return(scaled_df)
}


gtf_chm13<-"CHM13_hg38chrY_ecoliK12.gff3"
txdb_chm13<-makeTxDbFromGFF(gtf_chm13,dbxrefTag = "gene_name")

#Load peaks for CellType, MF and SKNF1 samples :
colnamesPeaks <- c("chr","start","end","name","score","strand","fold_enrichment","-log10pvalue","-log10qvalue","peak.point.source")
control_peaks <- fread("/gathered_peaks/Control_peaks_all.narrowPeak", sep="\t", col.names = colnamesPeaks)
kd_peaks <- fread("/gathered_peaks/Knockdown_peaks_all.narrowPeak", sep="\t", col.names = colnamesPeaks)
  
scaledPeaks_CellType_ctrlsh <- overlapPeakAndGenomeCoords2(control_peaks,makeGRangesFromDataFrame(metagenes_chm13_long,keep.extra.columns = T), tot.bins = tot.bins)

#####
# Number of peaks bar plot

peaksCellType <- rbind(
control_peaks %>% as.data.frame() %>% mutate(cond="Control"),
kd_peaks %>% as.data.frame() %>% mutate(cond="Knockdown"),
)

peaksCellType_subtel <- rbind(
  join_overlap_inner(control_peaks %>% 
                       mutate(arm=if_else(end > 100000,"Q","P"), strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
                       makeGRangesFromDataFrame(.,keep.extra.columns = T),chrEndCoords) %>% 
    as.data.frame() %>% filter((arm=="P" & strand=="-") | (arm=="Q" & strand=="+")) %>% mutate(cond="Control"),
  join_overlap_inner(kd_peaks %>% 
                       mutate(arm=if_else(end > 100000,"Q","P"), strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
                       makeGRangesFromDataFrame(.,keep.extra.columns = T),chrEndCoords) %>% 
    as.data.frame() %>% filter((arm=="P" & strand=="-") | (arm=="Q" & strand=="+")) %>% mutate(cond="KD")
)

ord1 <- c("Control","KD")

peaksCellType %>% group_by(cond) %>% dplyr::count() %>% ggplot(aes(fct_relevel(cond,ord1),n, color=cond)) + geom_col(width = 0.85, alpha=0.1,fill="white", lwd=1) + theme_prism(axis_text_angle = 45) + scale_y_continuous(labels = addUnitstext, expand = expansion(mult = c(0,0.1)), breaks = seq(0,25000,by=5000), guide = "prism_offset_minor") + labs(y="Number of peaks",x="",title="CellType") + theme(aspect.ratio = 1.3, legend.position = "none") + geom_text(aes(label=n), vjust=-1.1, color="black", size=3,position = position_dodge2(width = 0.9)) -> pTotalPeaks_CellType

peaksCellType_subtel %>% group_by(cond) %>% dplyr::count() %>% ggplot(aes(fct_relevel(cond,ord1),n, color=cond)) + geom_col(width = 0.85, alpha=0.1,fill="white", lwd=1) + theme_prism(axis_text_angle = 45) + scale_y_continuous(labels = addUnitstext, expand = expansion(mult = c(0,0.1)), breaks = seq(0,30,by=10), guide = "prism_offset_minor") + labs(y="Number of peaks",x="",title="Peaks at 50kb Chromosome ends \n CellType") + theme(aspect.ratio = 1.3, legend.position = "none") + geom_text(aes(label=n), vjust=-1.1, color="black", size=3,position = position_dodge2(width = 0.9)) -> pSubtelPeaks_CellType


ggsave("totalPeaks_CellType.pdf", plot = pTotalPeaks_CellType)
ggsave("subtelPeaks_CellType_subtel.pdf", plot = pSubtelPeaks_CellType)


### Convert to GRanges including strand information
control_peaks <- control_peaks %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

kd_peaks <- kd_peaks %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

#####
#Perparing metagene coordinated from CHM13 genome annotation
genesGTF_chm13 <- genes(txdb_chm13, columns=c("GENEID")) # extracted gene coordinates include 5'UTR and 3'UTR
genesGTF_chm13 <- genesGTF_chm13 %>% filter(!grepl("(NC_045512v2|U00096.2|J|G|M|K|chrUn|random|chrM)",seqnames,perl = T))

# Get transcript coordinates and select longest isoforms:
txdf_chm13<-AnnotationDbi::select(txdb_chm13,keys(txdb_chm13,"GENEID"),"TXNAME","GENEID")


txdf_chm13 <- left_join(transcripts(txdb_chm13) %>% as.data.frame(), txdf_chm13 ,by=c("tx_name"="TXNAME")) %>% 
  group_by(GENEID) %>% top_n(1,wt = width) %>% filter(., rank(dplyr::desc(width), ties.method = "first")==1) %>% ungroup() %>% distinct()

#chrY gene names
txdf_chm13 <- txdf_chm13 %>% mutate(tx_name=if_else(is.na(tx_name),GENEID,tx_name))

#####
#Get 5'UTR and 3'UTR coordinates, 

threeUTR_chm13 <- threeUTRsByTranscript(txdb_chm13, use.names=TRUE) %>% as.data.frame()

threeUTR_chm13 <- inner_join(threeUTR_chm13,txdf_chm13, by=c("group_name"="tx_name")) %>%
  dplyr::rename(threeUTR_start=start.x,threeUTR_end=end.x, threeUTR_strand=strand.x,tx_start=start.y,tx_end=end.y,tx_strand=strand.y) %>% 
  filter((tx_strand=="-" & threeUTR_start==tx_start) | (tx_strand=="+" & threeUTR_end==tx_end))

threeUTR_chm13 <- threeUTR_chm13 %>% dplyr::rename(threeUTR_seqnames=seqnames.x,threeUTR_width=width.x,tx_width=width.y,tx_seqnames=seqnames.y) %>% 
  dplyr::select(-c("group","exon_name","exon_id","exon_rank","tx_id"))

fiveUTR_chm13 <- fiveUTRsByTranscript(txdb_chm13, use.names=TRUE) %>% as.data.frame()

fiveUTR_chm13 <- inner_join(fiveUTR_chm13,txdf_chm13, by=c("group_name"="tx_name")) %>% 
  dplyr::rename(fiveUTR_start=start.x,fiveUTR_end=end.x, fiveUTR_strand=strand.x,tx_start=start.y,tx_end=end.y,tx_strand=strand.y) %>% 
  filter((tx_strand=="-" & fiveUTR_end==tx_end) | (tx_strand=="+" & fiveUTR_start==tx_start))

fiveUTR_chm13 <- fiveUTR_chm13 %>% dplyr::select(seqnames.x,fiveUTR_start,fiveUTR_end, width.x,fiveUTR_strand,group_name,GENEID) %>% 
  dplyr::rename(fiveUTR_seqnames=seqnames.x, fiveUTR_width=width.x)

# Extract CDS coordinates from 5'UTR, 3'UTR and transcript information and scale everything to generate metagene coordinates:
metagenes_chm13 <- inner_join(threeUTR_chm13,fiveUTR_chm13,by=c("GENEID"="GENEID","group_name"="group_name"))

metagenes_chm13 <- metagenes_chm13 %>% 
  mutate(cds_seqnames=tx_seqnames,
         cds_start=if_else(tx_strand=="+",fiveUTR_end+1,threeUTR_end+1),
         cds_end=if_else(tx_strand=="+",threeUTR_start-1,fiveUTR_start-1),
         cds_width=cds_end-cds_start,
         cds_strand=tx_strand,
         upstream_seqnames=tx_seqnames,
         upstream_start=if_else(tx_strand=="+",fiveUTR_start-1000,threeUTR_start-1000),
         upstream_end=if_else(tx_strand=="+",fiveUTR_start-1,threeUTR_start-1),
         upstream_width=upstream_end-upstream_start,
         upstream_strand=tx_strand,
         downstream_seqnames=tx_seqnames,
         downstream_start=if_else(tx_strand=="+",threeUTR_end+1,fiveUTR_end+1),
         downstream_end=if_else(tx_strand=="+",threeUTR_end+1000,fiveUTR_end+1000),
         downstream_width=downstream_end-downstream_start,
         downstream_strand=tx_strand,
  )


## Calculate bin sizes relative to cds sizes based on the median feature lengths:



getFeature<-function(df,label){
  
  res<- df %>% dplyr::select(starts_with(label),group_name,GENEID) %>% 
    mutate(type=label)
  
  colnames(res) <- gsub(paste0(label,"_"),"",colnames(res))
  return(res)
  
}

#Construct metagene object in long format
metagenes_chm13_long <- rbind(
  getFeature(metagenes_chm13,"upstream"),
  getFeature(metagenes_chm13,"fiveUTR"),
  getFeature(metagenes_chm13,"cds"),
  getFeature(metagenes_chm13,"threeUTR"),
  getFeature(metagenes_chm13,"downstream")
)

metagenes_chm13_long<- metagenes_chm13_long %>% mutate(mg_start=start,mg_end=end,mg_width=width,mg_strand=strand)


tot.bins<-3000
#####

genesGTF_chm13 <- genesGTF_chm13 %>% as.data.frame() %>% 
  mutate(gene_start=start,
         gene_end=end, 
         gene_width=gene_end-gene_start,
         start=if_else(gene_start >= 1000,gene_start-0,as.numeric(start)), #Expand gene coordinates +1kb around for metagene analysis
         end=if_else(gene_end >= 1000,gene_end+0,as.numeric(end)),
         mg_start=start,
         mg_end=end
  ) %>% makeGRangesFromDataFrame(., keep.extra.columns = T)


#####
#### Translate peak coordinates to metagene coordinates
library(scales)
scaledPeaks_CellType_ctrlsh <- overlapPeakAndGenomeCoords2(control_peaks,makeGRangesFromDataFrame(metagenes_chm13_long,keep.extra.columns = T), tot.bins = tot.bins)
scaledPeaks_CellType_kd <- overlapPeakAndGenomeCoords2(kd_peaks,makeGRangesFromDataFrame(metagenes_chm13_long,keep.extra.columns = T), tot.bins = tot.bins)

scaledPeaks_CellTypesamples <- rbind(
  scaledPeaks_CellType_ctrlsh %>% mutate(cond="CellType Control"),
  scaledPeaks_CellType_kd %>% mutate(cond="CellType KD")
)

scaledPeaks_CellTypesamples$cond <- factor(scaledPeaks_CellTypesamples$cond,c("CellType Control","CellType KD"))



scaledPeaks_CellTypesamples %>% ggplot(aes(x=scaled_summit_bin, color=cond)) + ggmulti::geom_density_(aes(y=..density..),as.mix = T) + facet_wrap(~cond, ncol = 2) + labs(y="Count") +
  scale_x_continuous(breaks=c(1,1.5,2,2.5,3,4.5,6,6.5,7,7.5,8),labels=c("","1kb","","5'UTR","","CDS","","3'UTR","","1kb","")) + 
  scale_y_continuous(expand = expansion(mult=c(0.01,.1))) +
  theme_prism(base_size = 10,base_rect_size = 1,axis_text_angle = 45,border = T) + theme(legend.position = "none") + labs(x="", title = "Peak density ", subtitle="")

#____________________________________________________________________________________________________________________________________________________________________#
##### xPore T2T genome-wide predicted sites

xptable_genwide_t2t <- fread("/nanoporeData/xpore_diffmod/diffmod.table", sep = ",",header = T)

xporeT2T <- xptable_genwide_t2t %>% filter( diff_mod_rate_CTRL_vs_KD >  0, 
                                            pval_CTRL_vs_KD < 0.05,
                                            #`mod_rate_CTRL-REP1` > 0.5,
                                            grepl("..A..",kmer)) %>% dplyr::select(id,position,kmer)


# Intersect k-mer coordinates with gene coordinates
#%>% mutate(start=position-3,end=position+2, strand=".",width=end-start, peak.point.source=position)
xporeT2T.gr <- left_join(xporeT2T,
                         genesGTF_chm13 %>% as.data.frame() %>% rownames_to_column("gene"),
                         by=c("id"="gene")) %>% dplyr::select(seqnames,start,end,width,strand,position,id,kmer) %>% filter(!is.na(seqnames)) %>% mutate(start=position-3, end=position+2,summit=position) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)

#####
scaledPeaks_xPoreT2T <- join_overlap_inner_within_directed(
  xporeT2T.gr,
  makeGRangesFromDataFrame(metagenes_chm13_long, keep.extra.columns = T),
  maxgap = 0,
  minoverlap = 1
) %>% as.data.frame() %>%
  mutate(
    scaled_start = if_else(strand == "+",
                           ((start - mg_start) / (mg_end - mg_start)) * tot.bins,
                           (1 - (
                             (start - mg_start) / (mg_end - mg_start)
                           )) * tot.bins),
    scaled_end = if_else(strand == "+",
                         (1 - (mg_end - end) / (mg_end - mg_start)) * tot.bins,
                         (((mg_end - end) / (mg_end - mg_start)
                         )) * tot.bins),
    scaled_summit = if_else(strand == "+",
                            ((summit - mg_start) / (mg_end - mg_start)) *
                              tot.bins,
                            (1 - (
                              (summit - mg_start) / (mg_end - mg_start)
                            )) * tot.bins)
  )

scaledPeaks_xPoreT2T <- scaledPeaks_xPoreT2T %>%
  mutate(
    scaled_summit_bin= if_else(
      type=="upstream",
      scales::rescale(scaled_summit, to=c(1,2),from = c(0,tot.bins)),
      if_else(type=="fiveUTR",
              scales::rescale(scaled_summit, to=c(2,3),from = c(0,tot.bins)),
              if_else(type=="cds", 
                      scales::rescale(scaled_summit, to=c(3,6),from = c(0,tot.bins)),
                      if_else(type=="threeUTR",
                              scales::rescale(scaled_summit, to=c(6,7),from = c(0,tot.bins)),
                              if_else(type=="downstream",
                                      scales::rescale(scaled_summit, to=c(7,8),from = c(0,tot.bins)),
                                      0
                              ))))))
#####


#____________________________________________________________________________________________________________________________________________________________________#

countPeaksPerBin <- function(df){
  binwidth=0.1
  bins <- seq(0,8, by=binwidth)
 
  df %>% 
    mutate(bin=cut(scaled_summit_bin,breaks=bins, include.lowest=F)) %>% 
    dplyr::group_by(cond,bin) %>% 
    add_count(bin) %>% 
    dplyr::select(scaled_summit_bin,bin,countsPerBin=n,cond) %>% 
    mutate(bin=as.numeric(gsub("(\\(|,.*\\])","",bin,perl = T))) %>%
    mutate(countsPerBin=as.numeric(countsPerBin)) %>% 
    ungroup() %>%
    dplyr::select(bin,countsPerBin,cond) %>% 
    distinct() %>% 
    pivot_wider(names_from = cond,values_from = countsPerBin) %>% 
    replace(is.na(.),0) %>%
    pivot_longer(names_to = "cond",values_to = "countsPerBin", -c("bin")) %>% 
    ungroup() -> countsPerBin_df
  
  return(countsPerBin_df)
}

library(ggalt)
library(Rcpp)

addUnitstext <- function(n) {
  labels <- ifelse(n < 1000, n,  # less than thousands
                   ifelse(n < 1e6, paste0(round(n/1e3, digits = 2), 'k'),  # in thousands
                          ifelse(n < 1e9, paste0(round(n/1e6,digits = 2), 'M'),  # in millions
                                 ifelse(n < 1e12, paste0(round(n/1e9,digits = 2), 'B'), # in billions
                                        ifelse(n < 1e15, paste0(round(n/1e12,digits = 2), 'T'), # in trillions
                                               'too big!'
                                        )))))
  return(labels)
}

plotPeaksPerBin <- function(df){
  max <- max(df$countsPerBin)/10
  
  df %>% 
  ggplot(aes(bin,countsPerBin, color=cond)) + 
  #geom_area(outline.type = "full") +
  geom_line() +
  #labs(y="Density") +
  scale_x_continuous(breaks=c(1,1.5,2,2.5,3,4.5,6,6.5,7,7.5,8),labels=c("","1kb","","5'UTR","","CDS","","3'UTR","","1kb","")) +
  annotate("rect",color="white",fill="blue",xmin = c(2,3,6),xmax = c(3,6,7),ymin = -(max*0.35),ymax = -(max*0.1), alpha=0.7) +
  annotate("rect",color="white",fill="black",xmin = c(1,7),xmax = c(2,8),ymin = -(max*0.25),ymax = -(max*0.1), alpha=0.7) +
  scale_y_continuous(labels = addUnitstext,expand = expansion(mult=c(0,.1))) +
  theme_prism(base_size = 10,base_rect_size = 0.5,axis_text_angle = 0,border = T) +
  theme(legend.position = c(0.85,0.8), axis.ticks.x = element_blank(), aspect.ratio = 0.75, legend.key.size = unit(0.5,"cm"), axis.text.x = element_text(vjust = -1,hjust = 0.5)) + 
  labs(x="", title = "", subtitle="") + coord_cartesian(ylim = c(0,NA), clip = "off") -> p
  
  return(p)
}


# Peaks per bin with density layout

plotPeaksPerBinDensity <- function(df){
  binwidth <- 0.1

 dens_df <- df %>% group_by(cond) %>% 
   arrange(bin) %>% 
   dplyr::slice(rep(row_number(),countsPerBin)) %>%
   ungroup()

 dens_df %>% 
  ggplot(aes(bin, color=cond,fill=cond)) + 
  geom_density(aes(y= binwidth* after_stat(count)), alpha=0.2) +
   scale_x_continuous(breaks=c(1.5,2.5,4.5,6.5,7.5),labels=c("1kb","5'UTR","CDS","3'UTR","1kb")) +
   scale_y_continuous(labels = addUnitstext,expand = expansion(mult=c(0,.1)), guide = "prism_offset_minor") -> p
 
 max <- max(ggplot_build(p)$layout$panel_params[[1]]$y.range) #Calculate max y-axis value to add plot annotations
 min <- min(ggplot_build(p)$layout$panel_params[[1]]$y.range)
 
 p <- p + 
   annotate("rect",color="white",fill="blue",xmin = c(2,3,6),xmax = c(3,6,7),ymin = -(max*0.005),ymax = -(max*0.025), alpha=0.7) +#-(max*0.015)
   annotate("rect",color="white",fill="black",xmin = c(1,7),xmax = c(2,8),ymin = -(max*0.015),ymax = -(max*0.005), alpha=1) + 
   theme_prism(base_size = 10,base_rect_size = 0.5,axis_text_angle = 0,border = F) +
   theme(legend.position = c(0.85,0.8), aspect.ratio = 0.75, 
         legend.key.size = unit(0.5,"cm"), axis.text.x = element_text(vjust = -0.75,hjust = 0.5), 
         axis.ticks.x=element_blank()) + 
   labs(x="",y="Density peaks per bin", title = "", subtitle="") + coord_cartesian(ylim = c(0,NA), clip = "off") -> p
 
 
 return(p)
}


countsPerBin_peaks_u2os <- scaledPeaks_CellTypesamples %>% countPeaksPerBin(.)

countsPerBin_xPoreT2T <- scaledPeaks_xPoreT2T %>% mutate(cond="Control CellType") %>% countPeaksPerBin(.)


 
######
### Metagene Plot (Number of Peaks per bin)

countsPerBin_peaks_u2os %>% plotPeaksPerBinDensity(.) + labs(title = "CellType") -> pNumPeaks_CellType

countsPerBin_xPoreT2T %>% plotPeaksPerBinDensity(.)  + scale_fill_manual(values = "steelblue") + scale_color_manual(values = "steelblue") + labs(title = "xPore T2T \n Differentially modified sites genome-wide \n NNANN", y="k-mers per bin") -> pNumPeaks_xPoreT2T


rbind(
  countsPerBin_xPoreT2T  %>% mutate(cond="xPore"),
  countsPerBin_peaks_u2os %>% filter(grepl("Control",cond)) %>% mutate(cond="m6A-RIP")
    ) %>% plotPeaksPerBinDensity(.) -> pNumPeaks_RIPandxPore


pNumPeaksDensity_all <- pNumPeaks_CellType | pNumPeaks_MFtumors | pNumPeaks_SKNF1 | pNumPeaks_SKNBE2

ggsave("metagenePlotDensity_PeaksPerBin_CellType.pdf",plot = pNumPeaks_CellType, width = 8, height = 6)
ggsave("figure_compilation/metagenePlotDensity_m6aRIP_and_xPore_overlay.pdf",plot = pNumPeaks_RIPandxPore, width = 8, height = 6)


######## Overlap peaks m6A-RIP vs Nanopore within TES en sites
intersectTES_RIP_xpore <- join_overlap_inner(
scaledPeaks_CellTypesamples %>% filter(dplyr::between(scaled_summit_bin,5.75,6.25)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),
scaledPeaks_xPoreT2T %>% mutate(start=start-100,end=end+100) %>% filter(dplyr::between(scaled_summit_bin,5.75,6.25)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)
)

library(VennDiagram)
library(venneuler)
nTESpeaks_RIP <- scaledPeaks_CellTypesamples %>% filter(dplyr::between(scaled_summit_bin,5.75,6.25)) %>% select(name) %>% unlist() %>% length()
nTESpeaks_xPore <- scaledPeaks_xPoreT2T %>% mutate(start=start-100,end=end+100) %>% filter(dplyr::between(scaled_summit_bin,5.75,6.25)) %>% mutate(siteID=row_number()) %>% select(siteID) %>% unlist() %>% length()

nintersectTES_RIP_xpore  <- intersectTES_RIP_xpore %>% as.data.frame() %>% nrow() %>% unlist()
#nTESpeaks_RIP-nintersectTES_RIP_xpore = 3446, #nTESpeaks_xPore-nintersectTES_RIP_xpore, #nintersectTES_RIP_xpore
pdf("figure_compilation/vennDiagram_xPore_m6ARIP_aroundTES.pdf")

v <- venneuler(c(m6ARIP=3446, xPore=5730,"m6ARIP&xPore"=946))
v$labels <- c("m6A-RIP\n3446","xPore\n5730")
plot(v, 
  main="peaks near TES, m6A-RIP vs xPore",
  col = hue_pal()(2)
  )

text(0.5,0.5,"946")

dev.off()

####### Correlation m6A-RIP and xPore coordinates

#### Distance between overlapping sites, xPore vs m6A-RIP (extending xPore k-mer sites 150bp around)
overlap_RIP_xPpore <- join_overlap_inner(
xporeT2T.gr %>% as.data.frame() %>% mutate(start=start-150,end=end+150) %>% rename(summit_xPore=summit) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),
control_peaks %>% mutate(summit_rip=start+peak.point.source) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)
)

### Save xPore T2T sites to bed format:
# %>% mutate(start=start-150,end=end+150)
xporeT2T.gr %>% as.data.frame()  %>% rename(summit_xPore=summit) %>% dplyr::select(seqnames,start,end,id,summit=position,strand) %>% fwrite(., "figure_compilation/tables/coordinates_xPoreT2T_diffMondAbove0_pval0.05.bed", col.names = F, sep="\t")

control_peaks %>% mutate(summit_rip=start+peak.point.source) %>% dplyr::select(chr,start,end,name,summit=summit_rip,strand) %>% fwrite(., "figure_compilation/tables/coordinates_m6aRIP_controlsh_CellType.bed", col.names = F, sep="\t")


#### Calculate distance to nearest m6A-RIP peak from xPore predicted sites:
# bedtools closest -D b -a coordinates_xPoreT2T_diffMondAbove0_pval0.05_sorted.bed -b coordinates_m6aRIP_controlsh_CellType_sorted.bed > distanceToClosest_xPore_to_m6aRIP.txt

  geom_histogram(data = overlap_RIP_xPpore %>% as.data.frame() %>%
                   mutate(distance=summit_xPore-summit_rip), aes(distance),binwidth = binwidth, alpha=0.3) +
  coord_cartesian(xlim = c(-600,600)) +
  theme_prism(axis_text_angle = 45) +
  theme(aspect.ratio = 1.5) +
  scale_y_continuous(expand = expansion(c(0,NA)), guide = "prism_offset_minor") +
  scale_x_continuous(labels = paste0(seq(-1000,1000,by=100),""),breaks = seq(-1000,1000,by=100)) +
  labs(title=str_wrap("Distance from summit between overlapping m6A-RIP and xPore sites",35), y="Number of overlapping sites",x="Distance between summits (bp)")

  
  distances_xporeT2T_to_m6aRIP <- join_nearest(xporeT2T.gr,
control_peaks %>% mutate(summit_rip=start+peak.point.source) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),
) %>% as.data.frame() %>%  mutate(distance=summit_rip-summit) 
  
  distances_xporeT2T_to_m6aRIP_random <- cbind(
    distances_xporeT2T_to_m6aRIP %>% dplyr::select(matches("_xpore")),
    distances_xporeT2T_to_m6aRIP %>% dplyr::select(!matches("_xpore")) %>% slice(sample(1:n()))
  ) %>% mutate(distance=summit_rip-summit_xpore)
  
  
  binwidth=49
  bins <- seq(-1000,1000, by=binwidth)
  
  df_dist <- distances_xporeT2T_to_m6aRIP %>% 
    mutate(bin=cut(distance,breaks=bins, include.lowest=F)) %>%
    dplyr::group_by(bin) %>% 
    add_count(bin) %>%
    mutate(countsPerBin=n) %>% 
    mutate(bin=as.numeric(gsub("(\\(|,.*\\])","",bin,perl = T))) %>%
    mutate(countsPerBin=as.numeric(countsPerBin)) %>% 
    ungroup() %>% 
    arrange(bin) %>%
    dplyr::select(bin,countsPerBin) %>% 
    distinct() %>% 
    dplyr::slice(rep(row_number(),countsPerBin)) %>%
    ungroup()
  
    ggplot() +
    #geom_histogram(data = distances_xporeT2T_to_m6aRIP, aes(distance),binwidth = 30, alpha=0.9) +
    geom_density(data = distances_xporeT2T_to_m6aRIP, aes(distance, color=str_wrap("xPore distance to nearest m6A-RIP peaks",25)), alpha=0.9) +
    #geom_density(data= distances_xporeT2T_to_m6aRIP_random, aes(distance,color=str_wrap("xPore distance to random m6A-RIP peaks",25)), alpha=0.9) +
    theme_prism(axis_text_angle = 45,base_size = 10) +
      theme(aspect.ratio = 1, legend.position = c(0.75,0.75)) +
      scale_y_continuous(expand = expansion(c(0,NA)), guide = "prism_offset_minor") +
      scale_x_continuous(labels = paste0(seq(-1000,1000,by=100),""),breaks = seq(-1000,1000,by=100),  limits = c(-1000,1000), expand = expansion(mult = 0.01)) +
      scale_color_manual(values = c("black","grey60")) +
      labs(title=str_wrap("Density \n Distance from xPore sites to nearest m6A-RIP summit \n genome-wide",35), y="Density",x="Distance between summits (bp)") -> plot_distance_RIP_xPore
    
    ggsave("figure_compilation/density_distanceNearestSite_xPore_to_m6ARIP_genomewide.pdf", plot = plot_distance_RIP_xPore, height = 6, width = 8)
    
  
  
    #rbind(distances_xporeT2T_to_m6aRIP %>% mutate(cond="Normal"),distances_xporeT2T_to_m6aRIP_random %>% mutate(cond="Random")) %>% ggplot(aes(x="",distance, color=cond)) + geom_boxplot() + geom_violin()
  


  
#####--------------------------------------------------------------------------------------------------------------------------------  
binwidth=30
bins <- seq(-1000,1000, by=binwidth)

df <- overlap_RIP_xPpore %>% as.data.frame() %>%
  mutate(distance=summit_xPore-summit_rip) %>%
  mutate(bin=cut(distance,breaks=bins, include.lowest=F)) %>%
  dplyr::group_by(bin) %>% 
  add_count(bin) %>%
  mutate(countsPerBin=n) %>% 
  mutate(bin=as.numeric(gsub("(\\(|,.*\\])","",bin,perl = T))) %>%
  mutate(countsPerBin=as.numeric(countsPerBin)) %>% 
  ungroup() %>% 
  arrange(bin) %>%
  dplyr::select(bin,countsPerBin) %>% 
  distinct() %>% 
  dplyr::slice(rep(row_number(),countsPerBin)) %>%
  ungroup()
  
nrows <- overlap_RIP_xPpore %>% as.data.frame() %>% nrow()  

overlap_RIP_xPpore_random <- cbind(
  overlap_RIP_xPpore %>% as.data.frame() %>% dplyr::select(seqnames,start,end,strand,position,id,summit_xPore),
  distances_xporeT2T_to_m6aRIP %>% as.data.frame() %>% dplyr::select(name_rip,summit_rip) %>% slice(sample(1:nrows))
) %>% mutate(distance=summit_rip-summit_xPore)

  ggplot() + 
  geom_density(data=overlap_RIP_xPpore %>% as.data.frame() %>%
                 mutate(distance=summit_xPore-summit_rip), aes(x=distance, color=str_wrap("xPore distance to nearest m6A-RIP peaks",25)), alpha=0.2) +
  #geom_density(data=overlap_RIP_xPpore_random, aes(x=distance, color=str_wrap("xPore distance to random m6A-RIP peaks",25)), alpha=0.2) +
  theme_prism(axis_text_angle = 45, base_size = 10) +
  theme(aspect.ratio = 1, legend.position = c(0.75,0.75)) +
  scale_y_continuous(expand = expansion(c(0.1,NA)), guide = "prism_offset_minor") +
  scale_x_continuous(labels = paste0(seq(-1000,1000,by=100),""),breaks = seq(-1000,1000,by=100), limits = c(-1000,1000)) +
  scale_color_manual(values = c("black","grey60")) +
  labs(title=str_wrap("Density \n Distance from xPore sites to nearest m6A-RIP summit \n Overlapping m6A-RIP and xPore sites",35), y="Density",x="Distance between summits (bp)") -> plot_distance_RIP_xPore_overlap

  ggsave("figure_compilation/density_distance_overlappingSites_m6ARIP_xPore.pdf", plot = plot_distance_RIP_xPore_overlap, height = 6, width = 8)

  
##### Correlation number of sites per gene, xPore vs m6A-RIP:
  
  genes_df <- genesGTF_chm13 %>% as.data.frame() %>% rownames_to_column("gene") %>% dplyr::select(gene)
  
  peaksPerGene_rip <- join_overlap_inner(control_peaks %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), 
                     genesGTF_chm13 %>% as.data.frame() %>% rownames_to_column("gene") %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)) %>% as.data.frame() %>% dplyr::select(gene,name) %>% distinct() %>% group_by(gene) %>% count() %>% ungroup()
  
  peaksPerGene_xpore <- xporeT2T %>% group_by(id) %>% count() %>% ungroup()
  
   #Bind peaks per gene from xPore and m6A-RIP datasets
  peaksPerGene_xPoreAndRIP <- left_join(left_join(
    genes_df,peaksPerGene_rip,by="gene"),
    peaksPerGene_xpore, by=c("gene"="id")
  ) %>% replace(is.na(.),0)
  
  peaksPerGene_xPoreAndRIP <- peaksPerGene_xPoreAndRIP %>% rename(peaksPerGene_rip=n.x,peaksPerGene_xpore=n.y)
  
  
#### Correlation coordinates genome-wide
  genesAndSummits_rip <- join_overlap_inner(
  genesGTF_chm13 %>% as.data.frame() %>% rownames_to_column("gene") %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),
  control_peaks %>% makeGRangesFromDataFrame(.,keep.extra.columns = T) ) %>% 
    as.data.frame() %>% mutate(summit=start+peak.point.source) %>%
    dplyr::select(gene,summit) %>% distinct()
  
  genesAndSummits_xpore <- left_join(xporeT2T, 
              genesGTF_chm13 %>% as.data.frame() %>% rownames_to_column("gene"), by=c("id"="gene")) %>% mutate(summit=position, gene=id) %>%
    dplyr::select(gene,summit) %>% distinct()

  summits_xpore_rip <- left_join(
  left_join(genes_df,genesAndSummits_rip, by="gene") %>% rename(summit_rip=summit),
  genesAndSummits_xpore, by="gene") %>% rename(summit_xpore=summit) %>% replace(is.na(.),0) %>% filter(!(summit_rip==0 & summit_xpore==0))
  
  summits_xpore_rip <- left_join(genesGTF_chm13 %>% as.data.frame() %>% rownames_to_column("gene") %>% dplyr::select(gene,seqnames,start,end)  %>% distinct(),
            summits_xpore_rip, by="gene") %>% filter(!(is.na(summit_xpore) & is.na(summit_rip))) %>% 
    mutate(chr=gsub("chr","",seqnames)) %>% 
    mutate(chr=if_else(chr=="X",23,as.numeric(chr))) %>% 
    #mutate(summit_xpore=as.numeric(paste0(chr,".",summit_xpore)), 
    #       summit_rip=as.numeric(paste0(chr,".",summit_rip))) %>% 
    filter(!is.na(summit_rip),!is.na(summit_xpore))
  

  
  summits_xpore_rip %>% ggplot(aes( x=summit_rip, y=summit_xpore)) +
    geom_point(size=0.01, alpha=0.5) +
    stat_cor() +
    geom_smooth(method = "lm", se = F) +
    theme_prism(base_size = 10,base_line_size = 0.5) + 
    scale_y_continuous(guide = "prism_offset_minor") + 
    scale_x_continuous(guide = "prism_offset_minor") +
    theme(aspect.ratio = 1) +
    labs(title = str_wrap("Genome-wide correlation of site coordinates xPore vs m6A-RIP summits",30), y="Coordinates xPore",x="Coordinates m6A-RIP") ->p
  
  ggsave("figure_compilation/coordinates_correlation_genomewide_xPoreVSm6aRIP_summits.pdf", width = 4, height = 4, plot = p)
  
  summits_xpore_rip %>% dplyr::select(-c("gene","start","end","chr","seqnames")) %>% correlate()
  
### Venn diagram m6A-RIP vs xPore sites
  library(eulerr)
  
  sites_xPore <- xporeT2T %>% as.data.frame() %>% mutate(site=paste0(id,position)) %>% dplyr::select(site) %>% unlist() %>% length()
  sites_m6aRIP_ctrlsh  <- length(unique(peaksControl$name)) 
  overlap_xPore_m6aRIP <- overlap_RIP_xPpore %>% as.data.frame() %>% mutate(site=paste0(id,position)) %>% dplyr::select(site) %>% distinct() %>% unlist() %>% length()
  
  vd <- euler(c("xPore Nanopore"=sites_xPore,"m6A-RIP"=sites_m6aRIP_ctrlsh,"xPore Nanopore&m6A-RIP"=overlap_xPore_m6aRIP))
  
  pdf("figure_compilation/vennDiagram_overlap_xPore_m6ARIP_sites.pdf", width = 6, height = 6)
  
  plot(vd,quantities=T,fills=c("steelblue","red"),lwd=NA,cex=0.1, alpha=0.3)
  
  dev.off()
  

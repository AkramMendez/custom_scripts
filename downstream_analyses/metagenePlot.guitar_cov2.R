library(tximport)
library(tximeta)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Guitar)
library(rtracklayer)
library(plyranges)
library(GenomicFeatures)
library(AnnotationDbi)


# Custom Functions

library(ggprism)

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
      #join_overlap_inner_within_directed(df,genesGTF,maxgap = 0,minoverlap = 1)
  scaled_df <-  join_overlap_inner(df,genesGTF,minoverlap = 1) %>% #join_overlap_inner_directed(df,genesGTF) # Will take peak summits with a min overlap of 1bp, in the same strand, annotated peaks that are not in the same strand or are farther away (more than 1kb) from any gene coordinate, will be discarded,
    as.data.frame() %>% 
    group_by(name) %>% top_n(1,wt = tx_width) %>% filter(., rank(desc(tx_width), ties.method = "first")==1) %>% ungroup() %>% mutate(GENEID=unlist(GENEID)) %>% distinct() %>%  # If a peak overlaps with several annotated genes, select the largest one.
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



gtf_chlsab2_cov2<-"/reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_concat.gtf"

gtf_grch38_cov2<-"/reference_genomes/custom_reference_genomes/gencode.v36.chr_patch_hapl_scaff.basic.annotation.gtf"

txdb.filename="/reference_genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_concat_txdb.sqlite"

txdb_chlsab2_cov2<-makeTxDbFromGFF(gtf_chlsab2_cov2,format = "gtf")

txdb_grch38_cov2<-makeTxDbFromGFF(gtf_grch38_cov2,format = "gtf")

#---- 
#Load peaks for samples :
colnamesPeaks <- c("chr","start","end","name","score","strand","fold_enrichment","-log10pvalue","-log10qvalue","peak.point.source")
LC_mock_d4 <- fread("/peakCalling/peaks_mock_d4_all.narrowPeak", sep="\t", col.names = colnamesPeaks)
LC_mock_d7 <- fread("/peakCalling/peaks_mock_d7_all.narrowPeak", sep="\t", col.names = colnamesPeaks)
LC_cov2_d4 <- fread("/peakCalling/peaks_cov2_d4_all.narrowPeak", sep="\t", col.names = colnamesPeaks) 
LC_cov2_d7 <- fread("/peakCalling/peaks_cov2_d7_all.narrowPeak", sep="\t", col.names = colnamesPeaks) 
N27_mock <- fread("/peakCalling/peaks_N27_mock_all.narrowPeak", sep="\t", col.names = colnamesPeaks) 
N27_cov2 <- fread("/peakCalling/peaks_N27_cov2_all.narrowPeak", sep="\t", col.names = colnamesPeaks)

N27_mock %>% filter(name %in% N27_pass) %>% as.data.frame() %>% fwrite(.,"peaks_HNE_mock.narrowPeak", sep="\t",col.names = F)

N27_cov2 %>% filter(name %in% N27_pass) %>% as.data.frame() %>% fwrite(.,"peaks_HNE_cov2.narrowPeak", sep="\t",col.names = F)



#Comparing pileup of samples
N27_mock_xls <- readPeaks_xls("/peakCalling/peaks_N27_mock_gathered.xls") %>% mutate(cond="N27 Mock")
N27_cov2_xls <- readPeaks_xls("/peakCalling/peaks_N27_cov2_gathered.xls") %>% mutate(cond="N27 SARS-CoV-2")

HBE_mock_d4_xls <- readPeaks_xls("/peakCalling/peaks_LC_mock_d4_gathered.xls") %>% mutate(cond="HBE Mock day 4")
HBE_cov2_d4_xls <- readPeaks_xls("/peakCalling/peaks_LC_cov2_d4_gathered.xls") %>% mutate(cond="HBE SARS-CoV-2 day 4")

N27_xls <- rbind(N27_mock_xls %>% as.data.frame(),N27_cov2_xls %>% as.data.frame())

HBE_d4_xls <- rbind(HBE_mock_d4_xls %>% as.data.frame(), HBE_cov2_d4_xls %>% as.data.frame())

LC_mock_d4_xls <- readPeaks_xls("/peakCalling/peaks_LC_mock_d4_gathered.xls")

rbind(LC_mock_d4_xls %>% as.data.frame() %>% mutate(cond="mock d4"),N27_cov2_xls %>% as.data.frame() %>% mutate(cond="N27 cov2")) %>% as.data.frame() %>% ggplot(aes(cond,pileup)) + geom_boxplot()

N27_xls %>% ggplot(aes(pileup, fold_enrichment)) + geom_point() + facet_wrap(~cond)
HBE_d4_xls %>% ggplot(aes(pileup, fold_enrichment)) + geom_point() + facet_wrap(~cond)



Calu_DE_DMSO <- fread("/peakCalling/peaks_Calu_samples_noDwnsmp/gathered_peaks/peaks_Calu3_DE_DMSO_gathered.narrowPeak", sep="\t", col.names = colnamesPeaks) %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)

Calu_DE_Sel <- fread("/peakCalling/peaks_Calu_samples_noDwnsmp/gathered_peaks/peaks_Calu3_DE_Sel_gathered.narrowPeak", sep="\t", col.names = colnamesPeaks) %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)

Calu_nonInf <- fread("/peakCalling/peaks_Calu_samples_noDwnsmp/gathered_peaks/peaks_Calu3_nonInf_gathered.narrowPeak", sep="\t", col.names = colnamesPeaks) %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)

Calu_UK_DMSO <- fread("/peakCalling/peaks_Calu_samples_noDwnsmp/gathered_peaks/peaks_Calu3_UK_DMSO_gathered.narrowPeak", sep="\t", col.names = colnamesPeaks) %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)

Calu_UK_Sel <- fread("/peakCalling/peaks_Calu_samples_noDwnsmp/gathered_peaks/peaks_Calu3_UK_Sel_gathered.narrowPeak", sep="\t", col.names = colnamesPeaks) %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)

################ Calu3 ext75

# Calu_DE_DMSO_ext75 <- fread("/peakCalling/peaks_Calu_samples_noDwnsmp/gathered_peaks/ext75/peaks_Calu3_DE_DMSO_gathered_ext75.narrowPeak", sep="\t", col.names = colnamesPeaks) %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)
# 
# Calu_DE_Sel_ext75 <- fread("/peakCalling/peaks_Calu_samples_noDwnsmp/gathered_peaks/ext75/peaks_Calu3_DE_Sel_gathered_ext75.narrowPeak", sep="\t", col.names = colnamesPeaks) %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)
# 
# Calu_nonInf_ext75 <- fread("/peakCalling/peaks_Calu_samples_noDwnsmp/gathered_peaks/ext75/peaks_Calu3_nonInf_gathered_ext75.narrowPeak", sep="\t", col.names = colnamesPeaks) %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)
# 
# Calu_UK_DMSO_ext75 <- fread("/peakCalling/peaks_Calu_samples_noDwnsmp/gathered_peaks/ext75/peaks_Calu3_UK_DMSO_gathered_ext75.narrowPeak", sep="\t", col.names = colnamesPeaks) %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)
# 
# Calu_UK_Sel_ext75<- fread("/peakCalling/peaks_Calu_samples_noDwnsmp/gathered_peaks/ext75/peaks_Calu3_UK_Sel_gathered_ext75.narrowPeak", sep="\t", col.names = colnamesPeaks) %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)



LC_mock_d4 <- LC_mock_d4 %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

LC_mock_d7 <- LC_mock_d7 %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

LC_cov2_d4 <- LC_cov2_d4 %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

LC_cov2_d7 <- LC_cov2_d7 %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

N27_mock <- N27_mock %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

N27_cov2 <- N27_cov2 %>% mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

vero_peaks.gr <- vero_peaks %>% 
  mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

wu_peaks.gr <- wu_peaks %>% 
  mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

uk_peaks.gr <- uk_peaks %>% 
  mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

sa_peaks.gr <- sa_peaks %>% 
  mutate(strand=if_else(grepl("_fwd_",name),"+","-")) %>% 
  makeGRangesFromDataFrame(.,keep.extra.columns = T)



genesGTF_cov2 <- genes(txdb_chlsab2_cov2, columns=c("GENEID")) # extracted gene coordinates include 5'UTR and 3'UTR
genesGTF_cov2 <- genesGTF_cov2 %>% filter(!grepl("(NC_045512v2|U00096.2|chrUn|random|chrM)",seqnames,perl = T))

genesGTF_hs <- genes(txdb_grch38_cov2, columns=c("GENEID")) # extracted gene coordinates include 5'UTR and 3'UTR
genesGTF_hs <- genesGTF_hs %>% filter(!grepl("(NC_045512v2|U00096.2|J|G|M|K|chrUn|random|chrM)",seqnames,perl = T))

# Get transcript coordinates and select longest isoforms:
txdf_hs<-AnnotationDbi::select(txdb_grch38_cov2,keys(txdb_grch38_cov2,"GENEID"),"TXNAME","GENEID")


txdf_hs <- left_join(transcripts(txdb_grch38_cov2) %>% as.data.frame(), txdf_hs ,by=c("tx_name"="TXNAME")) #%>% 
  #group_by(GENEID) %>% top_n(1,wt = width) %>% filter(., rank(desc(width), ties.method = "first")==1) %>% ungroup() %>% distinct()  


txdf_chlsab <-AnnotationDbi::select(txdb_chlsab2_cov2,keys(txdb_chlsab2_cov2,"GENEID"),"TXNAME","GENEID")

txdf_chlsab <- left_join(transcripts(txdb_chlsab2_cov2) %>% as.data.frame(), txdf_chlsab ,by=c("tx_name"="TXNAME")) #%>% 
  #group_by(GENEID) %>% top_n(1,wt = width) %>% filter(., rank(desc(width), ties.method = "first")==1) %>% ungroup() %>% distinct()  

#Get 5'UTR and 3'UTR coordinates, 

#### Homo sampiens


threeUTR_hs <- threeUTRsByTranscript(txdb_grch38_cov2, use.names=TRUE) %>% as.data.frame()

threeUTR_hs <- inner_join(threeUTR_hs,txdf_hs, by=c("group_name"="tx_name")) %>%
  rename(threeUTR_start=start.x,threeUTR_end=end.x, threeUTR_strand=strand.x,tx_start=start.y,tx_end=end.y,tx_strand=strand.y) %>% 
  filter((tx_strand=="-" & threeUTR_start==tx_start) | (tx_strand=="+" & threeUTR_end==tx_end))

threeUTR_hs <- threeUTR_hs %>% rename(threeUTR_seqnames=seqnames.x,threeUTR_width=width.x,tx_width=width.y,tx_seqnames=seqnames.y) %>% 
  dplyr::select(-c("group","exon_name","exon_id","exon_rank","tx_id"))

fiveUTR_hs <- fiveUTRsByTranscript(txdb_grch38_cov2, use.names=TRUE) %>% as.data.frame()
fiveUTR_hs <- inner_join(fiveUTR_hs,txdf_hs, by=c("group_name"="tx_name")) %>% 
  rename(fiveUTR_start=start.x,fiveUTR_end=end.x, fiveUTR_strand=strand.x,tx_start=start.y,tx_end=end.y,tx_strand=strand.y) %>% 
  filter((tx_strand=="-" & fiveUTR_end==tx_end) | (tx_strand=="+" & fiveUTR_start==tx_start))

fiveUTR_hs<- fiveUTR_hs %>% dplyr::select(seqnames.x,fiveUTR_start,fiveUTR_end, width.x,fiveUTR_strand,group_name,GENEID) %>% 
  rename(fiveUTR_seqnames=seqnames.x, fiveUTR_width=width.x)

##### Green monkey

threeUTR_chlsab <- threeUTRsByTranscript(txdb_chlsab2_cov2, use.names=TRUE) %>% as.data.frame()

threeUTR_chlsab <- inner_join(threeUTR_chlsab,txdf_chlsab, by=c("group_name"="tx_name")) %>%
  rename(threeUTR_start=start.x,threeUTR_end=end.x, threeUTR_strand=strand.x,tx_start=start.y,tx_end=end.y,tx_strand=strand.y) %>% 
  filter((tx_strand=="-" & threeUTR_start==tx_start) | (tx_strand=="+" & threeUTR_end==tx_end))

threeUTR_chlsab <- threeUTR_chlsab %>% rename(threeUTR_seqnames=seqnames.x,threeUTR_width=width.x,tx_width=width.y,tx_seqnames=seqnames.y) %>% 
  dplyr::select(-c("group","exon_name","exon_id","exon_rank","tx_id"))

fiveUTR_chlsab <- fiveUTRsByTranscript(txdb_chlsab2_cov2, use.names=TRUE) %>% as.data.frame()
fiveUTR_chlsab <- inner_join(fiveUTR_chlsab,txdf_chlsab, by=c("group_name"="tx_name")) %>% 
  rename(fiveUTR_start=start.x,fiveUTR_end=end.x, fiveUTR_strand=strand.x,tx_start=start.y,tx_end=end.y,tx_strand=strand.y) %>% 
  filter((tx_strand=="-" & fiveUTR_end==tx_end) | (tx_strand=="+" & fiveUTR_start==tx_start))

fiveUTR_chlsab<- fiveUTR_chlsab %>% dplyr::select(seqnames.x,fiveUTR_start,fiveUTR_end, width.x,fiveUTR_strand,group_name,GENEID) %>% 
  rename(fiveUTR_seqnames=seqnames.x, fiveUTR_width=width.x)


# Extract CDS coordinates from 5'UTR, 3'UTR and transcript information and scale everything to generate metagene coordinates:
metagenes_hs <- inner_join(threeUTR_hs,fiveUTR_hs,by=c("GENEID"="GENEID","group_name"="group_name"))


gtf_all_hs <- rtracklayer::import(gtf_grch38_cov2)

noAnn_genes_hs <- gtf_all_hs %>% as.data.frame() %>% 
  filter(!(transcript_id %in% unique(metagenes_hs$group_name)) & type=="transcript") %>%
  mutate(group_name=transcript_id,
         threeUTR_seqnames=seqnames,
         threeUTR_start=if_else(strand=="+",end,as.integer(start-1000)),
         threeUTR_end=if_else(strand=="+",as.integer(end+1000),start),
         threeUTR_width=abs(threeUTR_end-threeUTR_start),
         threeUTR_strand=strand,
         tx_seqnames=seqnames,
         tx_start=start,
         tx_end=end,
         tx_width=width,
         tx_strand=strand,
         GENEID=gene_id,
         fiveUTR_seqnames=seqnames,
         fiveUTR_start=if_else(strand=="+",as.integer(start-500),end),
         fiveUTR_end=if_else(strand=="+",start,as.integer(end+500)),
         fiveUTR_width=abs(fiveUTR_end-fiveUTR_start),
         fiveUTR_strand=strand
  ) %>% dplyr::select(colnames(metagenes_hs))

metagenes_hs <- rbind(metagenes_hs,noAnn_genes_hs)

metagenes_hs <- metagenes_hs %>% 
  mutate(cds_seqnames=tx_seqnames,
         cds_start=if_else(tx_strand=="+",fiveUTR_end+1,threeUTR_end+1),#fiveUTR_end+1,threeUTR_end+1),
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


##### Green monkey

gtf_all_cov2 <- rtracklayer::import(gtf_chlsab2_cov2)

# Genes without annotated 3UTR or 5UTR in the annotation GTF file:
# noAnnGenes_threeUTR_cov2 <- genesGTF_cov2 %>% as.data.frame() %>% filter(!(GENEID %in% unique(threeUTR_chlsab$GENEID))) %>% select(GENEID) %>% distinct() %>% unlist()
# noAnnGenes_fiveUTR_cov2 <- genesGTF_cov2 %>% as.data.frame() %>% filter(!(GENEID %in% unique(fiveUTR_chlsab$GENEID))) %>% select(GENEID) %>% distinct() %>% unlist()
# 
# names(noAnnGenes_threeUTR_cov2) <- NULL
# names(noAnnGenes_fiveUTR_cov2) <- NULL


metagenes_chlsab <- inner_join(threeUTR_chlsab,fiveUTR_chlsab,by=c("GENEID"="GENEID","group_name"="group_name")) #%>% filter(!grepl("(chrUn_|_random)",fiveUTR_seqnames) & !grepl("(chrUn_|_random)",threeUTR_seqnames))
  #inner_join(threeUTR_chlsab,fiveUTR_chlsab,by=c("GENEID"="GENEID","group_name"="group_name"))


#If a gene has no annotated 5UTR, include coordinates -500 bp from transcript start. mean(metagenes_chlsab$fiveUTR_width)=408
#If a gene has no annotated 5UTR, include coordinates +1000 bp from transcript end. mean(metagenes_chlsab$threeUTR_width)=1017
#"group_name"        "threeUTR_seqnames" "threeUTR_start"    "threeUTR_end"      "threeUTR_width"    "threeUTR_strand"   "tx_seqnames"       "tx_start"          "tx_end"           
#"tx_width"          "tx_strand"         "GENEID"            "fiveUTR_seqnames"  "fiveUTR_start"     "fiveUTR_end"       "fiveUTR_width"     "fiveUTR_strand"
noAnn_genes_cov2 <- gtf_all_cov2 %>% as.data.frame() %>% 
  filter(!(transcript_id %in% unique(metagenes_chlsab$group_name)) & type=="transcript") %>%
  mutate(group_name=transcript_id,
         threeUTR_seqnames=seqnames,
         threeUTR_start=if_else(strand=="+",end,as.integer(start-1000)),
         threeUTR_end=if_else(strand=="+",as.integer(end+1000),start),
         threeUTR_width=abs(threeUTR_end-threeUTR_start),
         threeUTR_strand=strand,
         tx_seqnames=seqnames,
         tx_start=start,
         tx_end=end,
         tx_width=width,
         tx_strand=strand,
         GENEID=gene_id,
         fiveUTR_seqnames=seqnames,
         fiveUTR_start=if_else(strand=="+",as.integer(start-500),end),
         fiveUTR_end=if_else(strand=="+",start,as.integer(end+500)),
         fiveUTR_width=abs(fiveUTR_end-fiveUTR_start),
         fiveUTR_strand=strand
         ) %>% dplyr::select(colnames(metagenes_chlsab))

metagenes_chlsab <- rbind(metagenes_chlsab,noAnn_genes_cov2)



# group_name threeUTR_seqnames threeUTR_start threeUTR_end threeUTR_width threeUTR_strand tx_seqnames tx_start tx_end tx_width tx_strand  GENEID fiveUTR_seqnames fiveUTR_start
# 1 XM_007993457.1              chr1           4620         4911            292               +        chr1     2628   4911     2284         + SCGB1C1             chr1          2628
# 2 XM_008004389.1              chr1           9933        10168            236               +        chr1     4916  10168     5253         +    ODF3             chr1          4916
# 3 XM_008011467.1              chr1           9933        10168            236               +        chr1     4916  10168     5253         +    ODF3             chr1          4916
# 4 XM_008018618.1              chr1          10167        10379            213               +        chr1     4916  10379     5464         +    ODF3             chr1          4916
# 5 XM_007977105.1              chr1          24638        25398            761               +        chr1    18597  25398     6802         +   RIC8A             chr1         18597
# 6 XM_007984396.1              chr1          24638        25398            761               +        chr1    18597  25398     6802         +   RIC8A             chr1         18597
# fiveUTR_end fiveUTR_width fiveUTR_strand
# 1        3657          1030              +
#   2        7312          2397              +

metagenes_chlsab <- metagenes_chlsab %>% 
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
  
  res<- df %>% dplyr::select(starts_with(label),group_name,GENEID,tx_width) %>% 
    mutate(type=label)
  
  colnames(res) <- gsub(paste0(label,"_"),"",colnames(res))
  return(res)
  
}

metagenes_hs_long <- rbind(
  getFeature(metagenes_hs,"upstream"),
  getFeature(metagenes_hs,"fiveUTR"),
  getFeature(metagenes_hs,"cds"),
  getFeature(metagenes_hs,"threeUTR"),
  getFeature(metagenes_hs,"downstream")
)

metagenes_hs_long<- metagenes_hs_long %>% mutate(mg_start=start,mg_end=end,mg_width=width,mg_strand=strand)

#### Green monkey
metagenes_chlsab_long <- rbind(
  getFeature(metagenes_chlsab,"upstream"),
  getFeature(metagenes_chlsab,"fiveUTR"),
  getFeature(metagenes_chlsab,"cds"),
  getFeature(metagenes_chlsab,"threeUTR"),
  getFeature(metagenes_chlsab,"downstream")
)

metagenes_chlsab_long<- metagenes_chlsab_long %>% mutate(mg_start=start,mg_end=end,mg_width=width,mg_strand=strand)





  
  tot.bins<-3000
#####
  #test 1
genesGTF_cov2 <- genesGTF_cov2 %>% as.data.frame() %>% 
  mutate(gene_start=start,
         gene_end=end, 
         gene_width=gene_end-gene_start,
         start=if_else(gene_start >= 1000,gene_start-0,as.numeric(start)), #Expand gene coordinates +1kb around for metagene analysis
         end=if_else(gene_end >= 1000,gene_end+0,as.numeric(end)),
         mg_start=start,
         mg_end=end
         ) %>% makeGRangesFromDataFrame(., keep.extra.columns = T)


genesGTF_hs <- genesGTF_hs %>% as.data.frame() %>% 
  mutate(gene_start=start,
         gene_end=end, 
         gene_width=gene_end-gene_start,
         start=if_else(gene_start >= 1000,gene_start-0,as.numeric(start)), #Expand gene coordinates +1kb around for metagene analysis
         end=if_else(gene_end >= 1000,gene_end+0,as.numeric(end)),
         mg_start=start,
         mg_end=end
  ) %>% makeGRangesFromDataFrame(., keep.extra.columns = T)





scaledPeaks_vero <- overlapPeakAndGenomeCoords2(vero_peaks.gr,makeGRangesFromDataFrame(metagenes_chlsab_long,keep.extra.columns = T), tot.bins = tot.bins)
scaledPeaks_wu <- overlapPeakAndGenomeCoords2(wu_peaks.gr,makeGRangesFromDataFrame(metagenes_chlsab_long,keep.extra.columns = T), tot.bins = tot.bins) #overlapPeakAndGenomeCoords(wu_peaks.gr,genesGTF_cov2)
scaledPeaks_uk <- overlapPeakAndGenomeCoords2(uk_peaks.gr,makeGRangesFromDataFrame(metagenes_chlsab_long,keep.extra.columns = T), tot.bins = tot.bins) #overlapPeakAndGenomeCoords(uk_peaks.gr,genesGTF_cov2)
scaledPeaks_sa <- overlapPeakAndGenomeCoords2(sa_peaks.gr,makeGRangesFromDataFrame(metagenes_chlsab_long,keep.extra.columns = T), tot.bins = tot.bins) #overlapPeakAndGenomeCoords(sa_peaks.gr,genesGTF_cov2)




df<- rbind(scaledPeaks_vero %>% mutate(cond="Non-infected"),
      scaledPeaks_wu %>% mutate(cond="B.1"),
      scaledPeaks_uk %>% mutate(cond="B.1.1.7"),
      scaledPeaks_sa %>% mutate(cond="B.1.351")
      ) 

binwidth <- 60
df %>% ggplot(aes(x=scaled_summit, color=cond)) + geom_line(stat = "bin",binwidth=binwidth) + ylim(0,500)
df %>% ggplot(aes(x=scaled_summit, color=cond)) + geom_density(aes(y=..scaled..)) + ylim(0,500)

df %>% ggplot(aes(x=scaled_summit, color=cond)) + ggmulti::geom_density_(aes(y=binwidth*..count..),bw=binwidth) + 
  geom_histogram(stat = "bin", binwidth = binwidth, alpha=0.1) + facet_wrap(~cond) + labs(y="Count")

df %>% ggplot(aes(x=scaled_summit, color=cond)) + ggmulti::geom_density_(aes(y=..density..),as.mix = T, bw=binwidth) + facet_wrap(~cond) + labs(y="Density")
######
#scaledPeaks_LCmockd7 <- overlapPeakAndGenomeCoords2(LC_mock_d7,genesGTF_hs, tot.bins = tot.bins)
scaledPeaks_LCmockd4 <- overlapPeakAndGenomeCoords2(LC_mock_d4,makeGRangesFromDataFrame(metagenes_hs_long %>% mutate(idx=row_number()),keep.extra.columns = T), tot.bins = tot.bins)
scaledPeaks_LCmockd7 <- overlapPeakAndGenomeCoords2(LC_mock_d7,makeGRangesFromDataFrame(metagenes_hs_long %>% mutate(idx=row_number()),keep.extra.columns = T), tot.bins = tot.bins)
scaledPeaks_LCcov2d4 <- overlapPeakAndGenomeCoords2(LC_cov2_d4,makeGRangesFromDataFrame(metagenes_hs_long %>% mutate(idx=row_number()),keep.extra.columns = T), tot.bins = tot.bins)
scaledPeaks_LCcov2d7 <- overlapPeakAndGenomeCoords2(LC_cov2_d7,makeGRangesFromDataFrame(metagenes_hs_long %>% mutate(idx=row_number()),keep.extra.columns = T), tot.bins = tot.bins)
scaledPeaks_N27mock <- overlapPeakAndGenomeCoords2(N27_mock,makeGRangesFromDataFrame(metagenes_hs_long %>% mutate(idx=row_number()),keep.extra.columns = T), tot.bins = tot.bins)
scaledPeaks_N27cov2 <- overlapPeakAndGenomeCoords2(N27_cov2,makeGRangesFromDataFrame(metagenes_hs_long %>% mutate(idx=row_number()),keep.extra.columns = T), tot.bins = tot.bins)

# scaledPeaks_LCmockd4_ext75 <- overlapPeakAndGenomeCoords2(LC_mock_d4_ext75,makeGRangesFromDataFrame(metagenes_hs_long,keep.extra.columns = T), tot.bins = tot.bins)
# scaledPeaks_LCmockd7_ext75 <- overlapPeakAndGenomeCoords2(LC_mock_d7_ext75,makeGRangesFromDataFrame(metagenes_hs_long,keep.extra.columns = T), tot.bins = tot.bins)
# scaledPeaks_LCcov2d4_ext75 <- overlapPeakAndGenomeCoords2(LC_cov2_d4_ext75,makeGRangesFromDataFrame(metagenes_hs_long,keep.extra.columns = T), tot.bins = tot.bins)
# scaledPeaks_LCcov2d7_ext75 <- overlapPeakAndGenomeCoords2(LC_cov2_d7_ext75,makeGRangesFromDataFrame(metagenes_hs_long,keep.extra.columns = T), tot.bins = tot.bins)
# scaledPeaks_N27mock_ext75 <- overlapPeakAndGenomeCoords2(N27_mock_ext75,makeGRangesFromDataFrame(metagenes_hs_long,keep.extra.columns = T), tot.bins = tot.bins)
# scaledPeaks_N27cov2_ext75 <- overlapPeakAndGenomeCoords2(N27_cov2_ext75,makeGRangesFromDataFrame(metagenes_hs_long,keep.extra.columns = T), tot.bins = tot.bins)



peaks_samples <- rbind(
scaledPeaks_LCmockd4 %>% mutate(cond="mock d4"),
scaledPeaks_LCmockd7 %>% mutate(cond="mock d7"),
scaledPeaks_LCcov2d4 %>% mutate(cond="cov2 d4"),
scaledPeaks_LCcov2d7%>% mutate(cond="cov2 d7"),
scaledPeaks_N27mock %>% mutate(cond="N27 mock"),
scaledPeaks_N27cov2 %>% mutate(cond="N27 cov2")
)

peaks_samples_ext75 <- rbind(
  scaledPeaks_LCmockd4_ext75 %>% mutate(cond="mock d4"),
  scaledPeaks_LCmockd7_ext75 %>% mutate(cond="mock d7"),
  scaledPeaks_LCcov2d4_ext75 %>% mutate(cond="cov2 d4"),
  scaledPeaks_LCcov2d7_ext75 %>% mutate(cond="cov2 d7"),
  scaledPeaks_N27mock_ext75 %>% mutate(cond="N27 mock"),
  scaledPeaks_N27cov2_ext75 %>% mutate(cond="N27 cov2")
)


df_hs %>% ggplot(aes(x=scaled_summit, color=cond)) + ggmulti::geom_density_(aes(y=..density..),as.mix = T, bw=binwidth) + facet_wrap(~cond) + labs(y="Density")
df_hs %>% ggplot(aes(x=scaled_summit_bin, color=cond)) + ggmulti::geom_density_(aes(y=..density..),as.mix = T) + facet_wrap(~cond) + labs(y="Density")

df_hs %>% ggplot(aes(x=scaled_summit_bin, color=cond)) + ggmulti::geom_density_(aes(y=..density..),as.mix = T) + facet_wrap(~cond) + labs(y="Density") +
  scale_x_continuous(breaks=c(1,1.5,2,2.5,3,4.5,6,6.5,7,7.5,8),labels=c("","1kb","","5'UTR","","CDS","","3'UTR","","1kb",""))

df_hs %>% ggplot(aes(scaled_summit_bin, color=cond)) + ggmulti::geom_density_(aes(y=..count..)) + 
  geom_histogram(stat = "bin", binwidth = binwidth, alpha=0.1) + 
  scale_x_continuous(breaks=c(0,3000),labels=c("TSS","TES")) +
  facet_wrap(~cond) + labs(y="Count")
################################

peaks_samples$cond <- factor(peaks_samples$cond,c("mock d4","cov2 d4","mock d7","cov2 d7","N27 mock","N27 cov2"))
peaks_samples_ext75$cond <- factor(peaks_samples_ext75$cond,c("mock d4","cov2 d4","mock d7","cov2 d7","N27 mock","N27 cov2"))

peaks_samples %>% ggplot(aes(x=scaled_summit_bin, color=cond)) + ggmulti::geom_density_(aes(y=..density..),as.mix = T) +
  facet_wrap(~cond, ncol = 2) + 
  labs(y="Count") +
  scale_x_continuous(breaks=c(1,1.5,2,2.5,3,4.5,6,6.5,7,7.5,8),labels=c("","1kb","","5'UTR","","CDS","","3'UTR","","1kb","")) + 
  scale_y_continuous(expand = expansion(mult=c(0.01,.1))) +
  theme_prism(base_size = 10,base_rect_size = 1,axis_text_angle = 45,border = T) + theme(legend.position = "none") + labs(x="", title = "Peak count ", subtitle="no downsamp, default ext")


#### Greeen monkey

scaledPeaks_vero1 <- overlapPeakAndGenomeCoords2(vero_peaks.gr,makeGRangesFromDataFrame(metagenes_chlsab_long %>% mutate(idx=row_number()),keep.extra.columns = T), tot.bins = tot.bins)
scaledPeaks_wu <- overlapPeakAndGenomeCoords2(wu_peaks.gr,makeGRangesFromDataFrame(metagenes_chlsab_long %>% mutate(idx=row_number()),keep.extra.columns = T), tot.bins = tot.bins)
scaledPeaks_uk <- overlapPeakAndGenomeCoords2(uk_peaks.gr,makeGRangesFromDataFrame(metagenes_chlsab_long %>% mutate(idx=row_number()),keep.extra.columns = T), tot.bins = tot.bins)
scaledPeaks_sa <- overlapPeakAndGenomeCoords2(sa_peaks.gr,makeGRangesFromDataFrame(metagenes_chlsab_long %>% mutate(idx=row_number()),keep.extra.columns = T), tot.bins = tot.bins)

#peaks_Vero_b1
scaledPeaks_vero <- rbind(
  scaledPeaks_vero1 %>% mutate(cond="Non-infected"),
  scaledPeaks_wu %>% mutate(cond="B.1"),
  scaledPeaks_uk %>% mutate(cond="B.1.1.7"),
  scaledPeaks_sa %>% mutate(cond="B.1.351")
)

#peaks_Vero_b1 
scaledPeaks_vero %>%  ggplot(aes(x=scaled_summit_bin, color=cond)) + ggmulti::geom_density_(aes(y=..count..),as.mix = T) + facet_wrap(~cond) + labs(y="Count") +
  scale_x_continuous(breaks=c(1,1.5,2,2.5,3,4.5,6,6.5,7,7.5,8),labels=c("","1kb","","5'UTR","","CDS","","3'UTR","","1kb","")) +
  coord_cartesian(xlim = c(1,8)) +
  scale_y_continuous(expand = expansion(mult=c(0.01,.1))) +
  theme_prism(base_size = 10,base_rect_size = 1,axis_text_angle = 45,border = T) + theme(legend.position = "none") + labs(x="", title = "Peak count Vero samples ", subtitle="no downsamp, default ext")


# scaledDensities_veroSamples %>%  ggplot(aes(x=x, color=cond)) + ggmulti::geom_density_(aes(y=..count..),as.mix = T) + facet_wrap(~cond) + labs(y="Count") +
#   scale_x_continuous(breaks=c(1,1.5,2,2.5,3,4.5,6,6.5,7,7.5,8),labels=c("","1kb","","5'UTR","","CDS","","3'UTR","","1kb","")) + 
#   scale_y_continuous(expand = expansion(mult=c(0.01,.1))) +
#   theme_prism(base_size = 10,base_rect_size = 1,axis_text_angle = 45,border = T) + theme(legend.position = "none") + labs(x="", title = "Peak count Vero samples ", subtitle="no downsamp, default ext")

###### Spike-in scaled density plots for HBE samples
sf_LCmockd4 <- 0.75
sf_LCmockd7 <- 0.68
sf_LCcov2d4 <- 1
sf_LCcov2d7 <- 0.69

sf_N27mock <- 1
sf_N27cov2 <- 0.92

total.rows<- 512*4
#,n=total.rows
density_LCmockd4 <- density(scaledPeaks_LCmockd4$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_LCmockd4*y)
density_LCmockd7 <- density(scaledPeaks_LCmockd7$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_LCmockd7*y)
density_LCcov2d4 <- density(scaledPeaks_LCcov2d4$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_LCcov2d4*y)
density_LCcov2d7 <- density(scaledPeaks_LCcov2d7$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_LCcov2d7*y)
density_N27mock <- density(scaledPeaks_N27mock %>% filter(name %in% N27_pass) %>% select(scaled_summit_bin) %>% unlist())[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_N27mock*y)
density_N27cov2 <- density(scaledPeaks_N27cov2 %>% filter(name %in% N27_pass) %>% select(scaled_summit_bin) %>% unlist())[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_N27cov2*y)

scaledDensities_samples <- rbind(
  density_LCmockd4 %>% mutate(cond="HBE mock day 4", label="HBE day 4"),
  density_LCmockd7 %>% mutate(cond="HBE mock day 7",label="HBE day 7"),
  density_LCcov2d4 %>% mutate(cond="HBE SARS-CoV-2 day 4",label="HBE day 4"),
  density_LCcov2d7 %>% mutate(cond="HBE SARS-CoV-2 day 7",label="HBE day 7"),
  density_N27mock %>% mutate(cond="N27 mock",label="N27"),
  density_N27cov2 %>% mutate(cond="N27 SARS-CoV-2",label="N27")
)

scaledDensities_samples$cond <- factor(scaledDensities_samples$cond, c("HBE mock day 4","HBE SARS-CoV-2 day 4", "HBE mock day 7", "HBE SARS-CoV-2 day 7", "N27 mock", "N27 SARS-CoV-2"))
# Human bronchial epithelial
library(unikn)
mycolors<-as.character(unlist(pal_unikn_pair[c(12,11,8,7)]))
#brewer.paired(8)[c(7,8,5,6,3,4)]

scaledDensities_samples %>%
  filter(!grepl("N27",cond)) %>% 
  ggplot(aes(x,y,color=cond)) + 
  geom_polygon(alpha=0.1, aes(fill=cond)) + 
  facet_wrap(~label, ncol = 3) + 
  labs(y="Density") +
  scale_x_continuous(breaks=c(1,1.5,2,2.5,3,4.5,6,6.5,7,7.5,8),labels=c("","1kb","","5'UTR","","CDS","","3'UTR","","1kb","")) + 
  scale_y_continuous(expand = expansion(mult=c(0.01,.1))) +
  annotate("rect",color="white",fill="blue",xmin = c(2,3,6),xmax = c(3,6,7),ymin = -0.04,ymax = -0.01, alpha=0.7) +
  annotate("rect",color="white",fill="black",xmin = c(1,7),xmax = c(2,8),ymin = -0.025,ymax = -0.04, alpha=0.7) +
  scale_color_manual(values = mycolorsLC) +
  scale_fill_manual(values = mycolorsLC) +
  coord_cartesian(ylim = c(0, 0.55),clip = "off") +
  theme_prism(base_size = 10,base_rect_size = 0.5,axis_text_angle = 0,border = T) + 
  theme(legend.position = "bottom", axis.ticks.x = element_blank(), aspect.ratio = 0.75, legend.key.size = unit(0.5,"cm"), 
        axis.text.x = element_text(vjust = -1,hjust = 0.5)) + 
  labs(x="", title = "") -> pPeakDensity_LC

ggsave("peakDensity_metagene_samples.pdf", height = 4, width = 8, plot = pPeakDensity_LC)




#### Scaled densities plot
#scaling factors:
sf_vero <- 1
sf_wu <-0.72
sf_uk <- 0.81
sf_sa <- 0.98


density_vero <- density(scaledPeaks_vero$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_vero*y)
density_wu <- density(scaledPeaks_wu$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_wu*y)
density_uk <- density(scaledPeaks_uk$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_uk*y)
density_sa <- density(scaledPeaks_sa$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_sa*y)

scaledDensities_veroSamples <- rbind(
  density_vero %>% mutate(cond="Non-infected"),
  density_wu %>% mutate(cond="B.1"),
  density_uk %>% mutate(cond="B.1.1.7"),
  density_sa %>% mutate(cond="B.1.351")
)


scaledDensities_veroSamples$cond <- factor(scaledDensities_veroSamples$cond,c("Non-infected","B.1","B.1.1.7","B.1.351"))

mycolors <- c("grey50", "#d7191c","#7ade95", "#fdae61")

#upstream c(1,2)
# "fiveUTR" c(2,3)
# "cds" c(3,6)
# "threeUTR" c(6,7)
# "downstream" c(7,8)
                                
#mycolors <- c("grey60",mycolors)
scaledDensities_veroSamples %>%
  ggplot(aes(x,scaled_density,color=cond)) + 
  geom_line(size=0.75) +
  #geom_polygon(alpha=0.1, size=0.65) +#,aes(fill=cond)) +
  labs(y="Density") +
  scale_x_continuous(breaks=c(1,1.5,2,2.5,3,4.5,6,6.5,7,7.5,8),labels=c("","1kb","","5'UTR","","CDS","","3'UTR","","1kb","")) +
  annotate("rect",color="white",fill="blue",xmin = c(2,3,6),xmax = c(3,6,7),ymin = -0.06,ymax = -0.01, alpha=0.7) +
  #annotate("rect",color="white",fill="black",xmin = c(1,7),xmax = c(2,8),ymin = -0.025,ymax = -0.05, alpha=0.7) +
  scale_y_continuous(expand = expansion(mult=c(0.01,.1))) +
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  coord_cartesian(ylim = c(0, 0.6),clip = "off") +
  theme_prism(base_size = 10,base_rect_size = 0.5,axis_text_angle = 0,border = F) +
  theme(legend.position = c(0.85,0.8), axis.ticks.x = element_blank(), aspect.ratio = 0.75, legend.key.size = unit(0.5,"cm"), axis.text.x = element_text(vjust = -1,hjust = 0.5)) + 
  labs(x="", title = "", subtitle="") -> pPeakDensity_vero
ggsave("peakDensity_metagene_vero.svg", height = 4, width = 6)

##### Scaling densities further for each group relative to the total data.#########################################################################################

#scaling factors:
sf_vero <- 1
sf_wu <-0.72
sf_uk <- 0.81
sf_sa <- 0.98

# Number of peaks overlapping metagene coordinates by condition
nDensity_vero <- nrow(scaledPeaks_vero1)
nDensity_wu <- nrow(scaledPeaks_wu)
nDensity_uk <- nrow(scaledPeaks_uk)
nDensity_sa <- nrow(scaledPeaks_sa)

# Total number of peaks overlapping metagene coordinates for all conditions
total_nDensity <- nDensity_vero + nDensity_wu + nDensity_uk + nDensity_sa

# Relative density proportion of each group:
proportion_vero <- nDensity_vero/total_nDensity
proportion_wu <- nDensity_wu/total_nDensity
proportion_uk <- nDensity_uk/total_nDensity
proportion_sa <- nDensity_sa/total_nDensity

densityProportional_vero <- density(scaledPeaks_vero$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_vero*y*proportion_vero)
densityProportional_wu <- density(scaledPeaks_wu$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_wu*y*proportion_wu)
densityProportional_uk <- density(scaledPeaks_uk$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_uk*y*proportion_uk)
densityProportional_sa <- density(scaledPeaks_sa$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_sa*y*proportion_sa)

scaledDensitiesProportional_veroSamples <- rbind(
  densityProportional_vero %>% mutate(cond="Non-infected"),
  densityProportional_wu %>% mutate(cond="B.1"),
  densityProportional_uk %>% mutate(cond="B.1.1.7"),
  densityProportional_sa %>% mutate(cond="B.1.351")
)


scaledDensitiesProportional_veroSamples$cond <- factor(scaledDensitiesProportional_veroSamples$cond,c("Non-infected","B.1","B.1.1.7","B.1.351"))

scaledDensitiesProportional_veroSamples %>%
  ggplot(aes(x,scaled_density,color=cond)) + 
  #geom_line(size=0.75) +
  geom_polygon(alpha=0.1, size=0.65,aes(fill=cond)) +
  labs(y="Density") +
  scale_x_continuous(breaks=c(1,1.5,2,2.5,3,4.5,6,6.5,7,7.5,8),labels=c("","1kb","","5'UTR","","CDS","","3'UTR","","1kb","")) +
  annotate("rect",color="white",fill="blue",xmin = c(2,3,6),xmax = c(3,6,7),ymin = -0.06,ymax = -0.01, alpha=0.7) +
  annotate("rect",color="white",fill="black",xmin = c(1,7),xmax = c(2,8),ymin = -0.025,ymax = -0.05, alpha=0.7) +
  scale_y_continuous(expand = expansion(mult=c(0,.1))) +
  scale_color_manual(values = mycolors) +
  scale_fill_manual(values = mycolors) +
  coord_cartesian(ylim = c(0, 0.6),clip = "off") +
  theme_prism(base_size = 10,base_rect_size = 0.5,axis_text_angle = 0,border = T) +
  theme(legend.position = c(0.85,0.8), axis.ticks.x = element_blank(), aspect.ratio = 0.75, legend.key.size = unit(0.5,"cm"), axis.text.x = element_text(vjust = -1,hjust = 0.5)) + 
  labs(x="", title = "", subtitle="") #-> pPeakDensityProportional_vero

ggsave("peakDensityProportionalToAllSamples_metagene_vero.pdf", height = 4, width = 6, plot = pPeakDensityProportional_vero)

################################################################################################################################################################################################################################################################################################################################################################
############################ Scaled density plots for samples ###############################################################################################################

sf_LCmockd4 <- 0.75
sf_LCmockd7 <- 0.68
sf_LCcov2d4 <- 1
sf_LCcov2d7 <- 0.69

#sf_N27mock <- 1
#sf_N27cov2 <- 0.92

# Number of peaks overlapping metagene coordinates by condition
nDensity_LCmockd4 <- nrow(scaledPeaks_LCmockd4)
nDensity_LCmockd7 <- nrow(scaledPeaks_LCmockd7)
nDensity_LCcov2d4 <- nrow(scaledPeaks_LCcov2d4)
nDensity_LCcov2d7 <- nrow(scaledPeaks_LCcov2d7)

# Total number of peaks overlapping metagene coordinates for all conditions
total_nDensity_<- nDensity_LCmockd4 + nDensity_LCmockd7 + nDensity_LCcov2d4 + nDensity_LCcov2d7

# Relative density proportion of each group:
proportion_LCmockd4 <- nDensity_LCmockd4/total_nDensity_LC
proportion_LCmockd7 <- nDensity_LCmockd7/total_nDensity_LC
proportion_LCcov2d4 <- nDensity_LCcov2d4/total_nDensity_LC
proportion_LCcov2d7 <- nDensity_LCcov2d7/total_nDensity_LC


densityProportional_LCmockd4 <- density(scaledPeaks_LCmockd4$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_LCmockd4*proportion_LCmockd4*y)
densityProportional_LCmockd7 <- density(scaledPeaks_LCmockd7$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_LCmockd7*proportion_LCmockd7*y)
densityProportional_LCcov2d4 <- density(scaledPeaks_LCcov2d4$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_LCcov2d4*proportion_LCcov2d4*y)
densityProportional_LCcov2d7 <- density(scaledPeaks_LCcov2d7$scaled_summit_bin)[c("x","y")] %>% as.data.frame() %>% mutate(scaled_density=sf_LCcov2d7*proportion_LCcov2d7*y)

binwidth=0.1
bins <- seq(0,8, by=binwidth)

makeBins <- function(df, bins){
  cond.label <- unique(df$cond)
  
  binned.df <- df  %>%
    mutate(bin=cut(scaled_summit_bin,breaks=bins, include.lowest=F)) %>% 
    dplyr::group_by(cond,bin) %>% 
    add_count(bin) %>% 
    dplyr::select(scaled_summit_bin,bin,countsPerBin=n,cond) %>% 
    mutate(bin=as.numeric(gsub("(\\(|,.*\\])","",bin,perl = T))) %>%
    mutate(countsPerBin=as.numeric(countsPerBin)) %>% 
    ungroup()
  
  bins.df <- data.frame(x=bins) 
  #Intersect peaks data with a fixed bin range dataset, to compare the same bins for all samples, as some peaks will not be overlappin with some bins across this fixed range as they don't have any peak there, the values will be set to 0, meaning that the bin is empty.
  binned.df <- full_join(bins.df,binned.df,by=c("x"="bin")) %>% 
    mutate(scaled_summit_bin=replace_na(scaled_summit_bin,0),
           countsPerBin=replace_na(countsPerBin,0),
           cond=replace_na(cond,cond.label)) %>% 
    dplyr::rename(bin=x)
  
  return(binned.df)
}

#scaledPeaks_LCmockd4 %>% mutate(cond="HBE mock day 4") %>% makeBins(.,bins=bins),
scaledPeaks_<- rbind(
  scaledPeaks_LCmockd4 %>% mutate(cond="HBE mock day 4") ,
  scaledPeaks_LCmockd7 %>% mutate(cond="HBE mock day 7") ,
  scaledPeaks_LCcov2d4 %>% mutate(cond="HBE SARS-CoV-2 day 4") ,
  scaledPeaks_LCcov2d7 %>% mutate(cond="HBE SARS-CoV-2 day 7")
)

scaledPeaks_LC$cond <- factor(scaledPeaks_LC$cond,c("HBE mock day 4","HBE mock day 7","HBE SARS-CoV-2 day 4","HBE SARS-CoV-2 day 7"))


tmp <- scaledPeaks_%>%
  mutate(bin=cut(scaled_summit_bin,breaks=bins, include.lowest=F)) %>% 
  dplyr::group_by(cond,bin) %>% 
  add_count(bin) %>% 
  dplyr::select(scaled_summit_bin,bin,countsPerBin=n,cond) %>% 
  mutate(bin=as.numeric(gsub("(\\(|,.*\\])","",bin,perl = T))) %>%
  mutate(countsPerBin=as.numeric(countsPerBin)) %>% 
  ungroup()

#ggplot(aes(bin,scaled_density, color=cond)) + 
tmp %>% ggplot(aes(bin,countsPerBin, color=cond)) + 
  geom_line() +
  #geom_polygon() +
  #geom_smooth(se=FALSE) +
  #labs(y="Density") +
  scale_x_continuous(breaks=c(1,1.5,2,2.5,3,4.5,6,6.5,7,7.5,8),labels=c("","1kb","","5'UTR","","CDS","","3'UTR","","1kb","")) +
  annotate("rect",color="white",fill="blue",xmin = c(2,3,6),xmax = c(3,6,7),ymin = -0.06,ymax = -0.01, alpha=0.7) +
  annotate("rect",color="white",fill="black",xmin = c(1,7),xmax = c(2,8),ymin = -0.025,ymax = -0.05, alpha=0.7) +
  scale_y_continuous(expand = expansion(mult=c(0,.1))) +
  #scale_color_manual(values = mycolorsLC) +
  #scale_fill_manual(values = mycolorsLC) +
  #coord_cartesian(ylim = c(0, 0.6),clip = "off") +
  theme_prism(base_size = 10,base_rect_size = 0.5,axis_text_angle = 0,border = T) +
  theme(legend.position = c(0.85,0.8), axis.ticks.x = element_blank(), aspect.ratio = 0.75, legend.key.size = unit(0.5,"cm"), axis.text.x = element_text(vjust = -1,hjust = 0.5)) + 
  labs(x="", title = "", subtitle="") 


scaledDensitiesProportional_samples <- rbind(
  densityProportional_LCmockd4 %>% mutate(cond="HBE mock day 4", label="HBE day 4"),
  densityProportional_LCmockd7 %>% mutate(cond="HBE mock day 7",label="HBE day 7"),
  densityProportional_LCcov2d4 %>% mutate(cond="HBE SARS-CoV-2 day 4",label="HBE day 4"),
  densityProportional_LCcov2d7 %>% mutate(cond="HBE SARS-CoV-2 day 7",label="HBE day 7")
  #density_N27mock %>% mutate(cond="N27 mock",label="N27"),
  #density_N27cov2 %>% mutate(cond="N27 SARS-CoV-2",label="N27")
)

scaledDensitiesProportional_samples$cond <- factor(scaledDensitiesProportional_samples$cond, c("HBE mock day 4","HBE SARS-CoV-2 day 4", "HBE mock day 7", "HBE SARS-CoV-2 day 7"))

# Human bronchial epithelial
library(unikn)
mycolors<-as.character(unlist(pal_unikn_pair[c(12,11,8,7)]))

scaledDensitiesProportional_samples %>%
  filter(!grepl("N27",cond)) %>% 
  ggplot(aes(x,scaled_density,color=cond)) + 
  geom_polygon(alpha=0.1, aes(fill=cond)) + 
  #facet_wrap(~label, ncol = 3) + 
  labs(y="Density") +
  scale_x_continuous(breaks=c(1,1.5,2,2.5,3,4.5,6,6.5,7,7.5,8),labels=c("","1kb","","5'UTR","","CDS","","3'UTR","","1kb","")) + 
  scale_y_continuous(expand = expansion(mult=c(0.01,.1))) +
  annotate("rect",color="white",fill="blue",xmin = c(2,3,6),xmax = c(3,6,7),ymin = -0.04,ymax = -0.01, alpha=0.7) +
  annotate("rect",color="white",fill="black",xmin = c(1,7),xmax = c(2,8),ymin = -0.025,ymax = -0.04, alpha=0.7) +
  #scale_color_manual(values = mycolorsLC) +
  #scale_fill_manual(values = mycolorsLC) +
  coord_cartesian(ylim = c(0, 0.55),clip = "off") +
  theme_prism(base_size = 10,base_rect_size = 0.5,axis_text_angle = 0,border = T) + 
  theme(legend.position = "bottom", axis.ticks.x = element_blank(), aspect.ratio = 0.75, legend.key.size = unit(0.5,"cm"), 
        axis.text.x = element_text(vjust = -1,hjust = 0.5)) + 
  labs(x="", title = "") #-> pPeakDensityProportional_LC

ggsave("peakDensityProportionalToAllSamples_metagene_samples.pdf", height = 4, width = 8, plot = pPeakDensityProportional_LC)


library(ggmulti)
#https://great-northern-diver.github.io/ggmulti/index.html
#https://github.com/great-northern-diver/ggmulti
#The ggmulti density function allows to calculate the grouped density using the following parameter, check help (?geom_density_) for more details un its use:
# as.mix	
# Logical. Within each group, if TRUE, the sum of the density estimate area is mixed and scaled to maximum 1. The area of each subgroup (in general, within each group one color represents one subgroup) is proportional to the count; if FALSE the area of each subgroup is the same, with maximum 1. See details.

#"each density estimates is proportional to the overall count."

# There are four combinations of scale.y and as.mix.
# 
# scale.y = "group" and as.mix = FALSE
# The density estimate area of each subgroup (represented by each color) within the same group is the same.
# 
# scale.y = "group" and as.mix = TRUE
# The density estimate area of each subgroup (represented by each color) within the same group is proportional to its own counts.
# 
# scale.y = "data" and as.mix = FALSE
# The sum of density estimate area of all groups is scaled to maximum of 1. and the density area for each group is proportional to the its count. Within each group, the area of each subgroup is the same.
# 
# scale.y = "data" and as.mix = TRUE
# The sum of density estimate area of all groups is scaled to maximum of 1 and the area of each subgroup (represented by each color) is proportional to its own count

scaledDensities_veroSamples %>% ggplot(aes(scaled_density,color=cond)) + geom_density()
  #geom_density_(scale.y = "group",as.mix = TRUE)

#peaks_Vero_b1
scaledPeaks_vero$cond <- factor(scaledPeaks_vero$cond, c("Non-infected","B.1","B.1.1.7","B.1.351"))

#peaks_Vero_b1 %>% 
scaledPeaks_vero %>% 
  ggplot(aes(x=scaled_summit_bin, color=cond)) + 
  ggmulti::geom_density_(as.mix = T, scale.y = "data") +  #aes(y=..density..)
  #facet_wrap(~cond) + labs(y="Count") +
  scale_x_continuous(breaks=c(1,1.5,2,2.5,3,4.5,6,6.5,7,7.5,8),labels=c("","1kb","","5'UTR","","CDS","","3'UTR","","1kb","")) + 
  scale_y_continuous(expand = expansion(mult=c(0,.1))) +
  scale_color_manual(values = mycolors) +
  theme_prism(base_size = 10,base_rect_size = 1,axis_text_angle = 45,border = T) +
  
#+ labs(x="", title = "Peak count Vero samples ", subtitle="no downsamp, default ext")


###############################################################################################################################################################################
#                                                                                                                                  
#                                                                    #UPDATED METAGENE PLOT USING PEAKS PER BIN COUNTS
#
###############################################################################################################################################################################
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
  
  dens_df <- df %>% 
    filter(dplyr::between(bin,2,7)) %>% 
    group_by(cond) %>% 
    arrange(bin) %>% 
    dplyr::slice(rep(row_number(),countsPerBin)) %>%
    ungroup()
  
  dens_df %>% 
    ggplot(aes(bin, color=cond,fill=cond)) + 
    geom_density(aes(y= binwidth* after_stat(count)), alpha=0.05, lwd=0.75) +#bw=binwidth #bw=binwidth/2
    #geom_histogram(alpha=0.3, position = "identity",binwidth = binwidth) +
    scale_x_continuous(breaks=c(2.5,4.5,6.5),labels=c("5'UTR","CDS","3'UTR")) +#scale_x_continuous(breaks=c(1.5,2.5,4.5,6.5,7.5),labels=c("1kb","5'UTR","CDS","3'UTR","1kb")) +
    scale_y_continuous(labels = addUnitstext,expand = expansion(mult=c(0,.1)), guide = "prism_offset_minor") -> p
  
  max <- max(ggplot_build(p)$layout$panel_params[[1]]$y.range) #Calculate max y-axis value to add plot annotations
  min <- min(ggplot_build(p)$layout$panel_params[[1]]$y.range)
  
  p <- p + 
    annotate("rect",color="white",fill="blue",xmin = c(2,3,6),xmax = c(3,6,7),ymin = -(max*0.005),ymax = -(max*0.025), alpha=0.7) +#-(max*0.015)
    #annotate("rect",color="white",fill="black",xmin = c(1,7),xmax = c(2,8),ymin = -(max*0.015),ymax = -(max*0.005), alpha=1) + 
    theme_prism(base_size = 10,base_rect_size = 0.5,axis_text_angle = 0,border = F) +
    theme(legend.position = c(0.85,0.85), aspect.ratio =1, 
          legend.key.size = unit(0.5,"cm"), axis.text.x = element_text(vjust = -0.75,hjust = 0.5), 
          axis.ticks.x=element_blank()) + 
    labs(x="",y="Peaks per bin", title = "", subtitle="") + coord_cartesian(ylim = c(0,NA), clip = "off") -> p
  
  
  return(p)
}


mycolors <- c("grey50", "#d7191c","#7ade95", "#fdae61")
mycolors<-c("#f1bbc8","#84a1e0","#a4ddd9","#fe894a")

peakAnns %>% filter(X.log10.qvalue._macs2 > qval_cutoff) %>% group_by(condition) %>% count()

commonLost <- intersect(intersect(lost_after_wu$name,lost_after_uk$name),lost_after_sa$name)
commonLostGenes <- intersect(intersect(lostGenes_wu$gene,lostGenes_uk$gene),lostGenes_sa$gene)

# Peaks per bin Vero and HBE samples:
countsPerBin_peaks_vero <- scaledPeaks_vero %>%
  #filter(GENEID %in% commonLostGenes) %>%
  #filter(name %in% union(union(gained_after_wu$name,gained_after_uk$name),gained_after_sa$name)) %>% 
  #filter(name %in% c(retained_after_sa$name,retained_after_uk$name,retained_after_wu$name)) %>%
  filter(!(name %in% c(fiverUTRStartPeaks_veroNonInf$name,fiverUTRStartPeaks_wuNonInf$name,fiverUTRStartPeaks_ukNonInf$name,fiverUTRStartPeaks_saNonInf))) %>% countPeaksPerBin(.) %>% ## Calculate counts peaks per bin and normalize to spike-in , remove peaks 20bp from TSS
  mutate(countsPerBin=if_else(grepl("Non-infected",cond), countsPerBin * sf_vero, 
                              if_else(grepl("B.1",cond), countsPerBin * sf_wu,
                                      if_else(grepl("B.1.1.7",cond), countsPerBin * sf_uk,
                                              if_else(grepl("B.1.351",cond), countsPerBin * sf_sa, countsPerBin
                                      )))))
#%>% filter(X.log10.q.value. > qval_cutoff)
countsPerBin_peaks_vero_noSpike <- scaledPeaks_vero %>% 
  #filter(GENEID %in% commonLostGenes) %>%
  #filter(name %in% union(union(gained_after_wu$name,gained_after_uk$name),gained_after_sa$name)) %>% 
  filter(!(name %in% c(fiverUTRStartPeaks_veroNonInf$name,fiverUTRStartPeaks_wuNonInf$name,fiverUTRStartPeaks_ukNonInf$name,fiverUTRStartPeaks_saNonInf))) %>% 
  countPeaksPerBin(.)

countsPerBin_peaks_samples <- scaledPeaks_%>% filter(!(name %in% c(fiveUTRStartPeaks_LCmockd4$name,fiveUTRStartPeaks_LCmockd7$name,fiveUTRStartPeaks_LCcov2d4$name,fiveUTRStartPeaks_LCcov2d7$name) )) %>% countPeaksPerBin(.) %>%  ## Calculate counts peaks per bin and normalize to spike-in
  mutate(countsPerBin=if_else(grepl("mock day 4",cond), countsPerBin  * sf_LCmockd4,
                              if_else(grepl("mock day 7",cond), countsPerBin * sf_LCmockd7,
                                      if_else(grepl("SARS-CoV-2 day 4",cond), countsPerBin * sf_LCcov2d4,
                                              if_else(grepl("SARS-CoV-2 day 7",cond), countsPerBin * sf_LCcov2d7,
                                                      countsPerBin
                                      )))))

# Peaks per bin density metagene plots (Vero samples):
countsPerBin_peaks_vero$cond <- factor(countsPerBin_peaks_vero$cond,c("Non-infected","B.1","B.1.1.7","B.1.351"))
countsPerBin_peaks_vero_noSpike$cond <- factor(countsPerBin_peaks_vero_noSpike$cond,c("Non-infected","B.1","B.1.1.7","B.1.351"))

countsPerBin_peaks_vero %>% plotPeaksPerBinDensity(.) + scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors) + theme(legend.position = c(0.9,0.9)) + labs(title = "Peaks per bin \n Spike-in Normalized")-> peaksPerBinPlot_Vero

countsPerBin_peaks_vero_noSpike %>% plotPeaksPerBinDensity(.) + scale_fill_manual(values = mycolors) + scale_color_manual(values = mycolors) + theme(legend.position = c(0.9,0.9)) +  labs(title = "Peaks per bin counts")-> peaksPerBinPlot_Vero_noSpike

p <- peaksPerBinPlot_Vero | peaksPerBinPlot_Vero_noSpike

# Count how many peaks are not included in the metagene and obtain its feature type:
scaledPeaks_vero %>% filter(cond=="Non-infected") %>% select(name) %>% distinct() %>% unlist() -> inMetagene_vero #12609
scaledPeaks_vero %>% filter(cond=="B.1") %>% select(name) %>% distinct() %>% unlist() -> inMetagene_wu #4954
scaledPeaks_vero %>% filter(cond=="B.1.1.7") %>% select(name) %>% distinct() %>% unlist() -> inMetagene_uk #4241
scaledPeaks_vero %>% filter(cond=="B.1.351") %>% select(name) %>% distinct() %>% unlist() -> inMetagene_sa #5854

peakAnns %>% filter(!(name_macs2 %in% inMetagene_vero) & condition=="vero" & X.log10.qvalue._macs2 > qval_cutoff) %>% count(short_annotation) #2251 Intergenic
peakAnns %>% filter(!(name_macs2 %in% inMetagene_wu) & condition=="wu" & X.log10.qvalue._macs2 > qval_cutoff) %>% count(short_annotation) #497 Intergenic
peakAnns %>% filter(!(name_macs2 %in% inMetagene_uk) & condition=="uk" & X.log10.qvalue._macs2 > qval_cutoff) %>% count(short_annotation) #1040 Intergenic
peakAnns %>% filter(!(name_macs2 %in% inMetagene_sa) & condition=="sa" & X.log10.qvalue._macs2 > qval_cutoff) %>% count(short_annotation) #2036 Intergenic

# Gained peaks in metagene
intersect(gained_after_wu$name,inMetagene_wu) %>% length() #2921
intersect(gained_after_uk$name,inMetagene_uk) %>% length() #2583
intersect(gained_after_sa$name,inMetagene_sa) %>% length() #3950

# Retained peaks in metagene
vero_wu %>% filter(`-log10(q.value)_vero` > qval_cutoff & name_wu %in% inMetagene_wu) %>% select(name_vero) %>% distinct() %>% dim()
vero_wu %>% filter(`-log10(q.value)_vero` > qval_cutoff & name_wu %in% inMetagene_wu) %>% select(name_wu) %>% distinct() %>% dim()# 2033 wu peaks correspond to 2107 peaks retained in vero
vero_uk %>% filter(`-log10(q.value)_vero` > qval_cutoff & name_uk %in% inMetagene_uk) %>% select(name_vero) %>% distinct() %>% dim() 
vero_uk %>% filter(`-log10(q.value)_vero` > qval_cutoff & name_uk %in% inMetagene_uk) %>% select(name_uk) %>% distinct() %>% dim() #1658 uk peaks correspond to 1763 retained peaks in vero
vero_sa %>% filter(`-log10(q.value)_vero` > qval_cutoff & name_sa %in% inMetagene_sa) %>% select(name_vero) %>% distinct() %>% dim()
vero_sa %>% filter(`-log10(q.value)_vero` > qval_cutoff & name_sa %in% inMetagene_sa) %>% select(name_sa) %>% distinct() %>% dim() #1904 sa peaks correspond to 2087 retained peaks in vero


#five UTR peaks
# Retained non-infected Vero
peakAnns %>% filter( condition=="vero" & X.log10.qvalue._macs2 > qval_cutoff) %>% select(name_macs2,condition) %>% filter(name_macs2 %in% fiverUTRStartPeaks_veroNonInf$name) %>%
distinct() %>% dim() #vero non-infected retained fiveUTR 20bp

retained_after_wu %>% filter(name %in% fiverUTRStartPeaks_veroNonInf$name) %>% select(name) %>% distinct() %>% dim() #6 retained after B.1 in fiveUTR 20bp
retained_after_uk %>% filter(name %in% fiverUTRStartPeaks_veroNonInf$name) %>% select(name) %>% distinct() %>% dim() #6 retained after B.1.1.7 in fiveUTR 20bp
retained_after_sa %>% filter(name %in% fiverUTRStartPeaks_veroNonInf$name) %>% select(name) %>% distinct() %>% dim() #3 retained after B.1.351 in fiveUTR 20bp

ret_wu_all <- retained_after_wu %>% select(name) %>% distinct() %>% unlist() 
ret_uk_all <- retained_after_uk %>% select(name) %>% distinct() %>% unlist() 
ret_sa_all <- retained_after_sa %>% select(name) %>% distinct() %>% unlist() 

vero_wu %>% filter(name_vero %in% ret_wu_all) %>% select(name_wu) %>% distinct() %>% dim() #2242 in Vero map to 2157 in WU, retained
vero_uk %>% filter(name_vero %in% ret_uk_all) %>% select(name_uk) %>% distinct() %>% dim() #1983 in Vero map to 1858 in UK, retained
vero_sa %>% filter(name_vero %in% ret_sa_all) %>% select(name_sa) %>% distinct() %>% dim() #2398 in Vero map to 2195 in SA, retained

ret_wu <- retained_after_wu %>% filter(name %in% fiverUTRStartPeaks_veroNonInf$name) %>% select(name) %>% distinct() %>% unlist()
ret_uk <- retained_after_uk %>% filter(name %in% fiverUTRStartPeaks_veroNonInf$name) %>% select(name) %>% distinct() %>% unlist()
ret_sa <- retained_after_sa %>% filter(name %in% fiverUTRStartPeaks_veroNonInf$name) %>% select(name) %>% distinct() %>% unlist()

vero_wu %>% filter(name_vero %in% ret_wu) %>% select(name_wu) %>% distinct() %>% dim() # 6 5UTR peaks retained in Vero, map to 6 peaks in WU
vero_uk %>% filter(name_vero %in% ret_uk) %>% select(name_uk) %>% distinct() %>% dim() # 6 5UTR peaks retained in Vero, map to 8 peaks in UK
vero_sa %>% filter(name_vero %in% ret_sa) %>% select(name_sa) %>% distinct() %>% dim() # 3 5UTR peaks retained in Vero, map to 5 peaks in SA

# Lost peaks at fiveUTR 20bp
lost_after_wu %>% filter(name %in% fiverUTRStartPeaks_veroNonInf$name) %>% select(name) %>% distinct() %>% dim() # 16 lost after B.1 in fiveUTR 20bp
lost_after_uk %>% filter(name %in% fiverUTRStartPeaks_veroNonInf$name) %>% select(name) %>% distinct() %>% dim() # 16 lost after B.1.1.7 in fiveUTR 20bp
lost_after_sa %>% filter(name %in% fiverUTRStartPeaks_veroNonInf$name) %>% select(name) %>% distinct() %>% dim() # 19 lost after B.1.351 in fiveUTR 20bp

# Gained peaks at fiveUTR 20bp
gained_after_wu %>% filter(name %in% fiverUTRStartPeaks_wuNonInf$name) %>% select(name) %>% distinct() %>% dim() # 10 gained after B.1 in fiveUTR 20bp
gained_after_uk %>% filter(name %in% fiverUTRStartPeaks_ukNonInf$name) %>% select(name) %>% distinct() %>% dim() # 2 gained after B.1.1.7 in fiveUTR 20bp
gained_after_sa %>% filter(name %in% fiverUTRStartPeaks_saNonInf$name) %>% select(name) %>% distinct() %>% dim() # 9 gained after B.1.351 in fiveUTR 20bp

ggsave("peaksPerBin_densityPlot_Vero_spikeInNorm_20221012.pdf", plot = peaksPerBinPlot_Vero, width = 6, height = 5, dpi = 300)

ggsave("peaksPerBin_densityPlot_Vero_spikeInNorm_and_noSpike_20221012.pdf", plot = p, width = 10, height = 5, dpi = 300)
ggsave("peaksPerBin_densityPlot_Vero_spikeInNorm_and_noSpike_CommonLostGenes_20221012.pdf", plot = p, width = 10, height = 5, dpi = 300)

ggsave("peaksPerBin_densityPlot_Vero_spikeInNorm.pdf", plot = peaksPerBinPlot_Vero, width = 6, height = 5, dpi = 300)

#ggsave("peaksPerBin_densityPlot_Vero_commonLostGenes.pdf", plot = peaksPerBinPlot_Vero, width = 6, height = 5, dpi = 300)


peakAnns %>% filter(X.log10.qvalue._macs2 > qval_cutoff) %>% mutate(in_metagene=if_else(name_macs2 %in% scaledPeaks_vero$name, TRUE,FALSE)) %>% count(condition,short_annotation, in_metagene) %>% filter(in_metagene==FALSE)


# Peaks per bin density metagene plots (HBE samples):
countsPerBin_peaks_samples$cond <- factor(countsPerBin_peaks_samples$cond, c("HBE mock day 4","HBE SARS-CoV-2 day 4","HBE mock day 7","HBE SARS-CoV-2 day 7"))
countsPerBin_peaks_samples %>% mutate(day=if_else(grepl("4",cond),"day 4","day 7")) %>% plotPeaksPerBinDensity(.) + facet_wrap(~day) + scale_fill_manual(values = mycolorsLC) + scale_color_manual(values = mycolorsLC) + theme(legend.position = c(0.7,0.75)) -> peaksPerBinPlot_HBEsamples

ggsave("peaksPerBin_densityPlot_HBEday4_HBEday7_spikeInNorm_20221017.pdf", plot = peaksPerBinPlot_HBEsamples, width = 8, height = 6, dpi = 300)

######## 
#Peaks per bin N27 samples
# Make barplot number of peaks, pileup above 3
N27_pass <- N27_xls %>% as.data.frame() %>% filter(pileup > 3 ) %>% select(name) %>% unlist()

peaks_samples %>% filter(name %in% N27_pass) %>%
  filter(!(name %in% c(fiveUTRStartPeaks_N27mock,fiveUTRStartPeaks_N27cov2))) %>% 
  mutate(cond=if_else(grepl("cov2",cond),"SARS-CoV-2","Mock")) %>% 
  countPeaksPerBin(.) %>% ### Calculate counts per bins and normalize to spike-in
  mutate(countsPerBin=if_else(grepl("Mock",cond), countsPerBin * sf_N27mock, countsPerBin * sf_N27cov2 )) %>% 
  plotPeaksPerBinDensity(.) + 
  scale_color_manual(values = pals::tol(4)[c(2,4)]) +
  scale_fill_manual(values = pals::tol(4)[c(2,4)]) ->  peaksPerBinPlot_N27samples
  
  ggsave("peaksPerBin_densityPlot_N27samples_pileupAbove3_spikeInNorm_20221017.pdf", plot = peaksPerBinPlot_N27samples, width = 8, height = 6, dpi = 300)
  
  
  p <- peaksPerBinPlot_Vero + labs(title = "Vero") | peaksPerBinPlot_HBEsamples + labs(title = "HBE") | peaksPerBinPlot_N27samples + labs(title = "NSE")
  
  ggsave("peaksPerBin_densityPlot_Vero_HBE_NSE_spikeInNorm_20221017.pdf", plot=p, width = 24, height = 6, dpi = 300)
  
  
###################################################################################################################  
########                                            Relative density plots  #########################################
  ###################################################################################################################    
library(ggmulti)  

plotDensity1 <- function(df){
   df %>%
    ggplot(aes(x,scaled_density,color=cond)) + 
      geom_line() +
      scale_x_continuous(breaks=c(2.5,4.5,6.5),labels=c("5'UTR","CDS","3'UTR")) +#scale_x_continuous(breaks=c(1.5,2.5,4.5,6.5,7.5),labels=c("1kb","5'UTR","CDS","3'UTR","1kb")) +
      scale_y_continuous(labels = addUnitstext,expand = expansion(mult=c(0,.1)), guide = "prism_offset_minor") -> p
    
    max <- max(ggplot_build(p)$layout$panel_params[[1]]$y.range) #Calculate max y-axis value to add plot annotations
    min <- min(ggplot_build(p)$layout$panel_params[[1]]$y.range)
    
    p <- p + 
      annotate("rect",color="white",fill="blue",xmin = c(2,3,6),xmax = c(3,6,7),ymin = -(max*0.005),ymax = -(max*0.025), alpha=0.7) +#-(max*0.015)
      #annotate("rect",color="white",fill="black",xmin = c(1,7),xmax = c(2,8),ymin = -(max*0.015),ymax = -(max*0.005), alpha=1) + 
      theme_prism(base_size = 10,base_rect_size = 0.5,axis_text_angle = 0,border = F) +
      theme(legend.position = c(0.85,0.85), aspect.ratio =1, 
            legend.key.size = unit(0.5,"cm"), axis.text.x = element_text(vjust = -0.75,hjust = 0.5), 
            axis.ticks.x=element_blank()) + 
      labs(x="",y="Peak Density", title="Spike-in Normalized", subtitle="") + coord_cartesian(ylim = c(0,NA), clip = "off") -> p
    
    
    return(p)
  }  
  
  
plotRelativeDensity_spike <- function(df){
  
      df %>%
    filter(dplyr::between(x,2,7)) %>%
      ggplot(aes(x,spikeInNorm_density ,color=cond,fill=cond)) +
      geom_line(lwd=0.75) +
      scale_x_continuous(breaks=c(2.5,4.5,6.5),labels=c("5'UTR","CDS","3'UTR")) +#scale_x_continuous(breaks=c(1.5,2.5,4.5,6.5,7.5),labels=c("1kb","5'UTR","CDS","3'UTR","1kb")) +
      scale_y_continuous(labels = addUnitstext,expand = expansion(mult=c(0,.1)), guide = "prism_offset_minor") -> p
    
    max <- max(ggplot_build(p)$layout$panel_params[[1]]$y.range) #Calculate max y-axis value to add plot annotations
    min <- min(ggplot_build(p)$layout$panel_params[[1]]$y.range)
    
    p <- p + 
      annotate("rect",color="white",fill="blue",xmin = c(2,3,6),xmax = c(3,6,7),ymin = -(max*0.005),ymax = -(max*0.025), alpha=0.7) +#-(max*0.015)
      #annotate("rect",color="white",fill="black",xmin = c(1,7),xmax = c(2,8),ymin = -(max*0.015),ymax = -(max*0.005), alpha=1) + 
      theme_prism(base_size = 10,base_rect_size = 0.5,axis_text_angle = 0,border = F) +
      theme(legend.position = c(0.85,0.85), aspect.ratio =1, 
            legend.key.size = unit(0.5,"cm"), axis.text.x = element_text(vjust = -0.75,hjust = 0.5), 
            axis.ticks.x=element_blank()) + 
      labs(x="",y="Relative Peak Density", title="Spike-in Normalized", subtitle="") + coord_cartesian(ylim = c(0,NA), clip = "off") -> p
    
    
    return(p)
  }


plotRelativeDensity_Nospike <- function(df){
  
  df %>% 
    filter(dplyr::between(x,2,7)) %>%
    ggplot(aes(x,y ,color=cond,fill=cond)) +
    geom_line(lwd=0.75) +
    scale_x_continuous(breaks=c(2.5,4.5,6.5),labels=c("5'UTR","CDS","3'UTR")) +#scale_x_continuous(breaks=c(1.5,2.5,4.5,6.5,7.5),labels=c("1kb","5'UTR","CDS","3'UTR","1kb")) +
    scale_y_continuous(labels = addUnitstext,expand = expansion(mult=c(0,.1)), guide = "prism_offset_minor") -> p
  
  max <- max(ggplot_build(p)$layout$panel_params[[1]]$y.range) #Calculate max y-axis value to add plot annotations
  min <- min(ggplot_build(p)$layout$panel_params[[1]]$y.range)
  
  p <- p + 
    annotate("rect",color="white",fill="blue",xmin = c(2,3,6),xmax = c(3,6,7),ymin = -(max*0.005),ymax = -(max*0.025), alpha=0.7) +#-(max*0.015)
    #annotate("rect",color="white",fill="black",xmin = c(1,7),xmax = c(2,8),ymin = -(max*0.015),ymax = -(max*0.005), alpha=1) + 
    theme_prism(base_size = 10,base_rect_size = 0.5,axis_text_angle = 0,border = F) +
    theme(legend.position = c(0.85,0.85), aspect.ratio =1, 
          legend.key.size = unit(0.5,"cm"), axis.text.x = element_text(vjust = -0.75,hjust = 0.5), 
          axis.ticks.x=element_blank()) + 
    labs(x="",y="Relative Peak Density",title= "No Normalization", subtitle="") + coord_cartesian(ylim = c(0,NA), clip = "off") -> p
  
  
  return(p)
}

scaledPeaks_vero$cond <- factor(scaledPeaks_vero$cond, c("Non-infected","B.1","B.1.1.7","B.1.351"))
  
  relDensity_vero <- scaledPeaks_vero %>% ggplot(aes(x=scaled_summit_bin, color=cond)) + ggmulti::geom_density_(aes(y=..density..),as.mix = T,scale.y = "data") + xlim(c(1,8))
  
          #F8766D= B.1, #C77CFF = noninf, #7CAE00=b.1.17, #00BFC4= B.1.351
  
  relDensity_vero_spikeIn <- ggplot_build(relDensity_vero)$data[[1]]
  
  relDensity_vero_spikeIn <- relDensity_vero_spikeIn %>% mutate(cond=if_else(colour=="#F8766D","Non-infected",if_else(colour=="#7CAE00","B.1",if_else(colour=="#00BFC4","B.1.1.7",if_else(colour=="#C77CFF","B.1.351",colour))))) ##C77CFF
  
  
  relDensity_vero_spikeIn <- relDensity_vero_spikeIn %>% mutate(cond=if_else(colour=="#F8766D","Non-infected",if_else(colour=="#7CAE00","B.1",if_else(colour=="#00BFC4","B.1.1.7",if_else(colour=="#C77CFF","B.1.351",colour))))) %>% 
    mutate(spikeInNorm_density=if_else(
    cond=="Non-infected",y*sf_vero,
    if_else(cond=="B.1",y*sf_wu,if_else(cond=="B.1.1.7",y*sf_uk,if_else(cond=="B.1.351",y*sf_sa,y)))))
  
  relDensity_vero_spikeIn$cond <- factor(relDensity_vero_spikeIn$cond,c("Non-infected","B.1","B.1.1.7","B.1.351"))
  relDensity_vero_spikeIn$cond <- factor(relDensity_vero_spikeIn$cond,c("Non-infected","B.1","B.1.1.7","B.1.351"))
  
  prelDensity_vero_spikeIn <- relDensity_vero_spikeIn %>% plotRelativeDensity_spike(.) + scale_color_manual(values =  mycolors)
  prelDensity_vero_NospikeIn <- relDensity_vero_spikeIn %>% plotRelativeDensity_Nospike(.) + scale_color_manual(values = mycolors)
  
  p2 <- prelDensity_vero_spikeIn | prelDensity_vero_NospikeIn
  
  ggsave("relativeDensity_spikeInNorm_NoSpike_Vero_20221012.pdf", width = 10, height = 5, dpi=300)
  
  ###### First version of density plots without relative normalization: --------------------------------------------------------------------------------------------------------
  scaledDensities_veroSamples %>% plotDensity1(.) + scale_color_manual(values = mycolors) -> pDensity1_vero
  
  scaledDensities_samples %>% filter(grepl("day 4",cond)) %>% mutate(scaled_density=y) %>% plotDensity1(.) + scale_color_manual(values = mycolorsLC[c(1,2)]) + theme(legend.position = c(1,0.85)) -> pDensity1_HBEd4
  scaledDensities_samples %>% filter(grepl("day 7",cond)) %>% mutate(scaled_density=y) %>% plotDensity1(.) + scale_color_manual(values = mycolorsLC[c(3,4)]) + theme(legend.position = c(1,0.85)) -> pDensity1_HBEd7
  scaledDensities_samples %>% filter(grepl("N27",cond)) %>% mutate(scaled_density=y) %>% plotDensity1(.) + scale_color_manual(values = pals::tol(4)[c(2,4)]) + theme(legend.position = c(1,0.85)) -> pDensity1_N27
  #################################################--------------------------------------------------------------------------------------------------------
  ############## Relative density plot HBE and N27
  
  fiveUTRpeaks_samples <- unique(c(fiveUTRStartPeaks_LCmockd4$name,fiveUTRStartPeaks_LCmockd7$name,fiveUTRStartPeaks_LCcov2d4$name,fiveUTRStartPeaks_LCcov2d7$name,fiveUTRStartPeaks_N27cov2$name,fiveUTRStartPeaks_N27mock$name))
  
  relDensity_HBEd4 <- peaks_samples %>% filter(grepl("d4",cond) & !(name %in% fiveUTRpeaks_samples)) %>% ggplot(aes(x=scaled_summit_bin, color=cond)) + ggmulti::geom_density_(aes(y=..density..),as.mix = T,scale.y = "data")
  
  relDensity_HBEd7 <- peaks_samples %>% filter(grepl("d7",cond)  & !(name %in% fiveUTRpeaks_samples)) %>% ggplot(aes(x=scaled_summit_bin, color=cond)) + ggmulti::geom_density_(aes(y=..density..),as.mix = T,scale.y = "data")
  
  relDensity_N27 <- peaks_samples %>% filter(grepl("N27",cond) & name %in% N27_pass  & !(name %in% fiveUTRpeaks_samples)) %>% ggplot(aes(x=scaled_summit_bin, color=cond)) + ggmulti::geom_density_(aes(y=..density..),as.mix = T,scale.y = "data")
  
  relDensity_HBEd4_spikeIn <- ggplot_build(relDensity_HBEd4)$data[[1]]
  relDensity_HBEd7_spikeIn <- ggplot_build(relDensity_HBEd7)$data[[1]]
  relDensity_N27_spikeIn <- ggplot_build(relDensity_N27)$data[[1]]
  
  
  #
  relDensity_HBEd4_spikeIn <- relDensity_HBEd4_spikeIn %>% mutate(cond=if_else(colour=="#00BFC4","SARS-CoV-2",if_else(colour=="#F8766D","Mock",colour))) %>% 
    mutate(spikeInNorm_density=if_else( cond=="Mock",y*sf_LCmockd4, if_else(cond=="SARS-CoV-2",y*sf_LCcov2d4,y)))
  
  relDensity_HBEd7_spikeIn <- relDensity_HBEd7_spikeIn %>% mutate(cond=if_else(colour=="#00BFC4","SARS-CoV-2",if_else(colour=="#F8766D","Mock",colour))) %>% 
    mutate(spikeInNorm_density=if_else( cond=="Mock",y*sf_LCmockd7, if_else(cond=="SARS-CoV-2",y*sf_LCcov2d7,y)))
  
  relDensity_N27_spikeIn <- relDensity_N27_spikeIn %>% mutate(cond=if_else(colour=="#00BFC4","SARS-CoV-2",if_else(colour=="#F8766D","Mock",colour)))  %>% 
    mutate(spikeInNorm_density=if_else( cond=="Mock",y*sf_N27mock, if_else(cond=="SARS-CoV-2",y*sf_N27cov2,y)))
  
  
  relDensity_HBEd4_spikeIn$cond <- factor(relDensity_HBEd4_spikeIn$cond,c("Mock","SARS-CoV-2"))
  relDensity_HBEd7_spikeIn$cond <- factor(relDensity_HBEd7_spikeIn$cond,c("Mock","SARS-CoV-2"))
  relDensity_N27_spikeIn$cond <- factor(relDensity_N27_spikeIn$cond,c("Mock","SARS-CoV-2"))
  
  
  prelDensity_HBEd4_spikeIn <- relDensity_HBEd4_spikeIn %>% plotRelativeDensity_spike(.) + scale_color_manual(values = mycolorsLC[c(1,2)]) + labs(subtitle = "HBE day 4")
  prelDensity_HBEd4_NospikeIn <- relDensity_HBEd4_spikeIn %>% plotRelativeDensity_Nospike(.) + scale_color_manual(values = mycolorsLC[c(1,2)])  + labs(subtitle = "HBE day 4")
  
  prelDensity_HBEd7_spikeIn <- relDensity_HBEd7_spikeIn %>% plotRelativeDensity_spike(.) + scale_color_manual(values = mycolorsLC[c(3,4)])  + labs(subtitle = "HBE day 7")
  prelDensity_HBEd7_NospikeIn <- relDensity_HBEd7_spikeIn %>% plotRelativeDensity_Nospike(.) + scale_color_manual(values = mycolorsLC[c(3,4)])  + labs(subtitle = "HBE day 7")
  
  prelDensity_N27_spikeIn <- relDensity_N27_spikeIn %>% plotRelativeDensity_spike(.)  + scale_color_manual(values = pals::tol(4)[c(2,4)]) + labs(subtitle = "NSE")
  prelDensity_N27_NospikeIn <- relDensity_N27_spikeIn %>% plotRelativeDensity_Nospike(.) + scale_color_manual(values = pals::tol(4)[c(2,4)]) + labs(subtitle = "NSE")
  
  prelDensity_HBE_N27 <- (prelDensity_HBEd4_spikeIn | prelDensity_HBEd4_NospikeIn) / ((prelDensity_HBEd7_spikeIn | prelDensity_HBEd7_NospikeIn)) / ((prelDensity_N27_spikeIn | prelDensity_N27_NospikeIn))
  
  ggsave("/plots/relativeDensity_HBE_N27_20221017.pdf", prelDensity_HBE_N27, width = 10, height = 16, dpi = 300)
  
  
  length(unique(intersecting_mock_cov2_LCd4$name_mock)) 
  length(unique(intersecting_mock_cov2_LCd4$name_mock)) 
  length(unique(gained_LCd4$name)) 
  
  length(unique(intersecting_mock_cov2_LCd7$name_cov2)) 
  length(unique(intersecting_mock_cov2_LCd7$name_mock)) 
  length(unique(gained_LCd7$name)) 
  
  # Retained peaks in the metagene, HBE day 4 and day 7, NSE
  peaks_samples %>% filter(name %in% unique(intersecting_mock_cov2_LCd4$name_mock) & cond=="mock d4") %>% select(name) %>% distinct() %>% dim() #191 retained in non-infected
  peaks_samples %>% filter(name %in% unique(intersecting_mock_cov2_LCd4$name_cov2)) %>% select(name) %>% distinct() %>% dim() # 178 retained from infected HBE d4
  
  peaks_samples %>% filter(name %in% unique(intersecting_mock_cov2_LCd7$name_mock) & cond=="mock d7") %>% select(name) %>% distinct() %>% dim() #669
  peaks_samples %>% filter(name %in% unique(intersecting_mock_cov2_LCd7$name_cov2)) %>% select(name) %>% distinct() %>% dim() # 579 retained from infected HBE d7
  #peaks_samples %>% filter(name %in% unique(intersecting_mock_cov2_LCd4$name_mock) & cond=="mock d7") %>% select(name) %>% distinct() %>% dim() #669
  
  # Gained peaks in the metagene, HBE day 4 and day 7, NSE
  peaks_samples %>% filter(name %in% unique(gained_LCd4$name) & cond=="cov2 d4") %>% select(name) %>% distinct() %>% dim() #367
  peaks_samples %>% filter(name %in% unique(gained_LCd7$name) & cond=="cov2 d7") %>% select(name) %>% distinct() %>% dim() #992
  
  retained_lost_gained_%>% filter(type=="Gained") %>% count(condition)
  # Retained, gained in metagene
  retained_lost_gained_%>% filter(name %in% unique(peaks_samples$name)) %>% select(condition,type,type2) %>% count(condition,type,type2)
  
  #Retained, lost, gained all HBE d4, d7
  retained_lost_gained_ %>% select(condition,type,type2) %>% count(condition,type,type2)
  
  retained_lost_gained_%>% count(condition,type2) 
  annPeaks_%>% count(sample)
  
  ret_from_mock_d7_5utr <- intersecting_mock_cov2_LCd7 %>% filter(name_mock %in% fiveUTRStartPeaks_LCmockd7$name) %>% select(name_cov2) %>% distinct() %>% unlist()
  
  library(hiAnnotator)
  peaks_table_<- getNearestFeature(annPeaks_%>% rename(seqnames=chr),gtf_all_hs,side = "either", colnam = "nearest_feature",feature.colnam = "gene_name", relativeTo = "subject") %>% as.data.frame()
  # Count HBE and NSE peaks not in the metagene
  peaks_table_%>% filter(peak_name %in% N27_pass | grepl("LC_",peak_name)) %>% filter(!(peak_name %in% c(peaks_samples$name))) %>% count(sample)

  # Count HBE and NSE peaks in the metagene
  peaks_table_%>% filter(peak_name %in% N27_pass | grepl("LC_",peak_name)) %>% filter(peak_name %in% c(peaks_samples$name)) %>% count(sample)
  
  peaks_table_%>% filter(peak_name %in% N27_pass | grepl("LC_",peak_name)) %>% filter(!(peak_name %in% c(peaks_samples$name))) %>% filter(!grepl("intergenic",short_annotation)) %>% summarise(range(nearest_featureDist))
  
  #ann_LCmockd7_fwd%>% filter(peak_name=="peaks_LC_mock_d7_fwd_peak_417")
  
  peaks_table_%>% filter(peak_name %in% N27_pass | grepl("LC_",peak_name)) %>% count(sample)
  
  peaks_table_%>% filter(peak_name %in% N27_pass | grepl("LC_",peak_name)) %>% filter(peak_name %in% c(peaks_samples$name)) %>% count(sample)
  
  ######### Plot density vs relative density Vero. HBE. N27
  
  pDensities_vero <- pDensity1_vero + labs(title="Spike-in Normalized \n Peak Density", y="Density") + theme(legend.position = "none") | prelDensity_vero_spikeIn + labs(title="Spike-in Normalized \n Relative Peak Density", y="Relative Density") + theme(legend.position = c(1,0.85))
  
  
  ggsave("/plots/peakDensity_vs_relativePeakDensity_vero.pdf", plot = pDensities_vero, width = 8, height = 4)
  
  pDensities_HBE_N27 <- (pDensity1_HBEd4  + labs(title="Spike-in Normalized \n Peak Density", subtitle = "HBE day 4", y="Density") + theme(legend.position = "none") | 
      prelDensity_HBEd4_spikeIn  + labs(title="Spike-in Normalized \n Relative Peak Density", subtitle = "HBE day 4", y="Relative Density")  + theme(legend.position = c(1,0.85))) / 
    (pDensity1_HBEd7  + labs(title="Spike-in Normalized \n Peak Density", subtitle = "HBE day 7", y="Density")  + theme(legend.position = "none") | 
       prelDensity_HBEd7_spikeIn  + labs(title="Spike-in Normalized \n Relative Peak Density",subtitle = "HBE day 7", y="Relative Density") + theme(legend.position = c(1,0.85))) / 
    (pDensity1_N27  + labs(title="Spike-in Normalized \n Peak Density", subtitle = "NSE", y="Density")  + theme(legend.position = "none") | 
       prelDensity_N27_spikeIn + labs(title="Spike-in Normalized \n Relative Peak Density",subtitle = "NSE", y="Relative Density") + theme(legend.position = c(1,0.85)) )
  
  
  ggsave("/plots/peakDensity_vs_relativePeakDensity_HBE_NSE.pdf", plot = pDensities_HBE_N27, width = 8, height = 10)

##########################
# 20 nt from CDS
##########################

  
overlapPeakAndGenomeCoords2(vero_peaks.gr,makeGRangesFromDataFrame(metagenes_chlsab_long %>% filter(type=="cds") %>% mutate(cds_start=start,cds_end=end,cds_strand=strand) %>% mutate(start=if_else(cds_strand=="+",cds_start,cds_end-20),end=if_else(cds_strand=="+",cds_start+20,cds_end)),keep.extra.columns = T), tot.bins = tot.bins)

# Vero
cdsStartPeaks_veroNonInf <- join_overlap_inner_within_directed(vero_peaks.gr %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), metagenes_chlsab_long %>% mutate(cds_seqnames=seqnames,cds_start=start,cds_end=end,cds_strand=strand,start=if_else(cds_strand=="+",cds_start,cds_end-20),end=if_else(cds_strand=="+",cds_start+20,cds_end)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()
  
cdsStartPeaks_B.1 <- join_overlap_inner_within_directed(wu_peaks.gr %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),metagenes_chlsab %>% mutate(start=if_else(cds_strand=="+",cds_start,cds_end-20),end=if_else(cds_strand=="+",cds_start+20,cds_end)) %>% makeGRangesFromDataFrame(.,seqnames.field = "cds_seqnames",keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()

cdsStartPeaks_B.1.1.7 <- join_overlap_inner_within_directed(uk_peaks.gr %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),metagenes_chlsab %>% mutate(start=if_else(cds_strand=="+",cds_start,cds_end-20),end=if_else(cds_strand=="+",cds_start+20,cds_end)) %>% makeGRangesFromDataFrame(.,seqnames.field = "cds_seqnames",keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()

cdsStartPeaks_B.1.351 <- join_overlap_inner_within_directed(sa_peaks.gr %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),metagenes_chlsab %>% mutate(start=if_else(cds_strand=="+",cds_start,cds_end-20),end=if_else(cds_strand=="+",cds_start+20,cds_end)) %>% makeGRangesFromDataFrame(.,seqnames.field = "cds_seqnames",keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()


##########################
# 20 nt from 5'UTR
##########################

fiverUTRStartPeaks_veroNonInf <- join_overlap_inner_within_directed(vero_peaks.gr %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), metagenes_chlsab_long %>% filter(type=="fiveUTR") %>% mutate(fiverUTR_seqnames=seqnames,fiverUTR_start=start,fiverUTR_end=end,fiverUTR_strand=strand,start=if_else(fiverUTR_strand=="+",fiverUTR_start,fiverUTR_end-20),end=if_else(fiverUTR_strand=="+",fiverUTR_start+20,fiverUTR_end)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()

fiverUTRStartPeaks_wuNonInf <- join_overlap_inner_within_directed(wu_peaks.gr %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), metagenes_chlsab_long %>% filter(type=="fiveUTR") %>% mutate(fiverUTR_seqnames=seqnames,fiverUTR_start=start,fiverUTR_end=end,fiverUTR_strand=strand,start=if_else(fiverUTR_strand=="+",fiverUTR_start,fiverUTR_end-20),end=if_else(fiverUTR_strand=="+",fiverUTR_start+20,fiverUTR_end)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()

fiverUTRStartPeaks_ukNonInf <- join_overlap_inner_within_directed(uk_peaks.gr %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), metagenes_chlsab_long %>% filter(type=="fiveUTR") %>% mutate(fiverUTR_seqnames=seqnames,fiverUTR_start=start,fiverUTR_end=end,fiverUTR_strand=strand,start=if_else(fiverUTR_strand=="+",fiverUTR_start,fiverUTR_end-20),end=if_else(fiverUTR_strand=="+",fiverUTR_start+20,fiverUTR_end)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()

fiverUTRStartPeaks_saNonInf <- join_overlap_inner_within_directed(sa_peaks.gr %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), metagenes_chlsab_long %>% filter(type=="fiveUTR") %>% mutate(fiverUTR_seqnames=seqnames,fiverUTR_start=start,fiverUTR_end=end,fiverUTR_strand=strand,start=if_else(fiverUTR_strand=="+",fiverUTR_start,fiverUTR_end-20),end=if_else(fiverUTR_strand=="+",fiverUTR_start+20,fiverUTR_end)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()


# HBE




fiveUTRStartPeaks_LCmockd4 <- join_overlap_inner_within_directed(LC_mock_d4 %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), metagenes_hs_long %>% filter(type=="fiveUTR") %>% mutate(fiveUTR_seqnames=seqnames,fiveUTR_start=start,fiveUTR_end=end,fiveUTR_strand=strand,start=if_else(fiveUTR_strand=="+",fiveUTR_start,fiveUTR_end-20),end=if_else(fiveUTR_strand=="+",fiveUTR_start+20,fiveUTR_end)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()

fiveUTRStartPeaks_LCmockd7 <- join_overlap_inner_within_directed(LC_mock_d7 %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),metagenes_hs_long %>% filter(type=="fiveUTR") %>% mutate(fiveUTR_seqnames=seqnames,fiveUTR_start=start,fiveUTR_end=end,fiveUTR_strand=strand,start=if_else(fiveUTR_strand=="+",fiveUTR_start,fiveUTR_end-20),end=if_else(fiveUTR_strand=="+",fiveUTR_start+20,fiveUTR_end)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()

fiveUTRStartPeaks_LCcov2d4 <- join_overlap_inner_within_directed(LC_cov2_d4 %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),metagenes_hs_long %>% filter(type=="fiveUTR") %>% mutate(fiveUTR_seqnames=seqnames,fiveUTR_start=start,fiveUTR_end=end,fiveUTR_strand=strand,start=if_else(fiveUTR_strand=="+",fiveUTR_start,fiveUTR_end-20),end=if_else(fiveUTR_strand=="+",fiveUTR_start+20,fiveUTR_end)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()

fiveUTRStartPeaks_LCcov2d7 <- join_overlap_inner_within_directed(LC_cov2_d7 %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),metagenes_hs_long %>% filter(type=="fiveUTR") %>% mutate(fiveUTR_seqnames=seqnames,fiveUTR_start=start,fiveUTR_end=end,fiveUTR_strand=strand,start=if_else(fiveUTR_strand=="+",fiveUTR_start,fiveUTR_end-20),end=if_else(fiveUTR_strand=="+",fiveUTR_start+20,fiveUTR_end)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()


fiveUTRStartPeaks_N27mock <- join_overlap_inner_within_directed(peaks_samples %>% filter(cond=="N27 mock" & name %in% N27_pass) %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),metagenes_hs_long %>% filter(type=="fiveUTR") %>% mutate(fiveUTR_seqnames=seqnames,fiveUTR_start=start,fiveUTR_end=end,fiveUTR_strand=strand,start=if_else(fiveUTR_strand=="+",fiveUTR_start,fiveUTR_end-20),end=if_else(fiveUTR_strand=="+",fiveUTR_start+20,fiveUTR_end)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()

fiveUTRStartPeaks_N27cov2 <- join_overlap_inner_within_directed(peaks_samples %>% filter(cond=="N27 cov2" & name %in% N27_pass) %>% mutate(startOriginal=start,start=(startOriginal+peak.point.source),end=(startOriginal+peak.point.source)+1) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),metagenes_hs_long %>% filter(type=="fiveUTR") %>% mutate(fiveUTR_seqnames=seqnames,fiveUTR_start=start,fiveUTR_end=end,fiveUTR_strand=strand,start=if_else(fiveUTR_strand=="+",fiveUTR_start,fiveUTR_end-20),end=if_else(fiveUTR_strand=="+",fiveUTR_start+20,fiveUTR_end)) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T), maxgap = 0, minoverlap = 1) %>% as.data.frame()


#### Type of five utr peaks:
    #Lost
    intersect(lost_after_wu$name,fiverUTRStartPeaks_veroNonInf$name) %>% unique() %>% length()
    intersect(lost_after_uk$name,fiverUTRStartPeaks_veroNonInf$name) %>% unique() %>% length()
    intersect(lost_after_sa$name,fiverUTRStartPeaks_veroNonInf$name) %>% unique() %>% length()
    #retained_lost_gained_%>% filter(name %in% fiveUTRStartPeaks_LCmockd4$name & type=="Lost" & type2=="Infected") %>% select(name) %>% distinct() %>% dim()
    #retained_lost_gained_%>% filter(name %in% fiveUTRStartPeaks_LCmockd7$name & type=="Lost" & type2=="Infected" ) %>% select(name) %>% distinct() %>% dim()
    
    lost_LCd4 %>% filter(name %in% fiveUTRStartPeaks_LCmockd4$name) %>% select(name) %>% distinct() %>% dim()
    lost_LCd7 %>% filter(name %in% fiveUTRStartPeaks_LCmockd7$name) %>% select(name) %>% distinct() %>% dim()
    lost_N27 %>% filter(name %in% fiveUTRStartPeaks_N27mock$name) %>% select(name) %>% distinct() %>% dim()
    

    #Gained
    intersect(gained_after_wu$name,fiverUTRStartPeaks_wuNonInf$name) %>% unique() %>% length()
    intersect(gained_after_uk$name,fiverUTRStartPeaks_ukNonInf$name) %>% unique() %>% length()
    intersect(gained_after_sa$name,fiverUTRStartPeaks_saNonInf$name) %>% unique() %>% length()
    
    #retained_lost_gained_%>% filter(name %in% fiveUTRStartPeaks_LCcov2d4$name & type=="Gained" & type2=="Infected")  %>% select(name) %>% distinct() %>% dim()
    #retained_lost_gained_%>% filter(name %in% fiveUTRStartPeaks_LCcov2d7$name & type=="Gained" & type2=="Infected" )  %>% select(name) %>% distinct() %>% dim()
    
    gained_LCd4 %>% filter(name %in% fiveUTRStartPeaks_LCcov2d4$name) %>% select(name) %>% distinct() %>% dim()
    gained_LCd7 %>% filter(name %in% fiveUTRStartPeaks_LCcov2d7$name) %>% select(name) %>% distinct() %>% dim()
    gained_N27 %>% filter(name %in% fiveUTRStartPeaks_N27cov2$name) %>% select(name) %>% distinct() %>% dim()

    # Retained
    vero_wu %>% filter(name_wu %in% fiverUTRStartPeaks_wuNonInf$name) %>% select(name_vero) %>% distinct() %>% dim()
    vero_uk %>% filter(name_uk %in% fiverUTRStartPeaks_ukNonInf$name) %>% select(name_vero) %>% distinct() %>% dim()
    vero_sa %>% filter(name_sa %in% fiverUTRStartPeaks_saNonInf$name) %>% select(name_vero) %>% distinct() %>% dim()
    
    retained_LCd4 %>% filter(name %in% fiveUTRStartPeaks_LCmockd4$name) %>% select(name) %>% distinct() %>% dim()
    retained_LCd7 %>% filter(name %in% fiveUTRStartPeaks_LCmockd7$name) %>% select(name) %>% distinct() %>% dim()
    retained_N27 %>% filter(name %in% fiveUTRStartPeaks_N27mock$name) %>% select(name) %>% distinct() %>% dim()
    
    
#################################################################################
    # Barplot number of peaks in metagene by feature
#################################################################################
    
    order.type <- c("5'UTR","CDS","stop codon","3'UTR")
    
    countsPerBin_peaks_samples %>%
      #filter(dplyr::between(bin,2,7)) %>% 
      mutate(type=if_else(dplyr::between(bin,0,2),"upstream",
                          if_else(dplyr::between(bin,2,3),"5'UTR",
                                  if_else(dplyr::between(bin,3,5.8),"CDS",
                                          if_else(dplyr::between(bin,5.9,6),"stop codon",
                                                  if_else(dplyr::between(bin,6.1,7),"3'UTR",
                                                          if_else(dplyr::between(bin,6.2,8),"downstream","NA"))))))) %>% 
      group_by(cond,type) %>% 
      mutate(n=round(sum(countsPerBin))) %>% 
      ungroup() %>%
      select(cond,type,n) %>% 
      distinct() %>% 
      filter(grepl("day 4",cond) &! grepl("stream",type) ) %>%
      mutate(cond=if_else(grepl("CoV-2",cond),"SARS-CoV-2","Mock")) %>% 
      ggplot(aes(fct_relevel(type,order.type),n,fill=cond)) + 
      geom_col(position = position_dodge2()) +
      theme_prism(base_size = 10, axis_text_angle = 45) +
      scale_fill_manual(values = mycolorsLC[1:2]) +
      scale_y_continuous(guide = "prism_offset", expand = expansion(mult = c(0,NA)), limits = c(0,3500)) +
      labs(y="Number of peaks", x="", title = "HBE day 4") +
      #geom_text(aes(label=n), position = position_dodge2(width = 0.9), vjust=-0.5, size=3) +
      theme(legend.position = c(0.75,0.9), aspect.ratio = 0.75) -> pNpeaks_HBEd4
    
    
    
    countsPerBin_peaks_samples %>%
      #filter(dplyr::between(bin,2,7)) %>% 
      mutate(type=if_else(dplyr::between(bin,0,2),"upstream",
        if_else(dplyr::between(bin,2,3),"5'UTR",
                          if_else(dplyr::between(bin,3,5.8),"CDS",
                                  if_else(dplyr::between(bin,5.9,6),"stop codon",
                                          if_else(dplyr::between(bin,6.1,7),"3'UTR",
                                                  if_else(dplyr::between(bin,6.2,8),"downstream","NA"))))))) %>% 
      group_by(cond,type) %>% 
      mutate(n=round(sum(countsPerBin))) %>% 
      ungroup() %>%
      select(cond,type,n) %>% 
      distinct() %>% 
      filter(grepl("day 7",cond) &! grepl("stream",type) ) %>%
      mutate(cond=if_else(grepl("CoV-2",cond),"SARS-CoV-2","Mock")) %>% 
      ggplot(aes(fct_relevel(type,order.type),n,fill=cond)) + 
      geom_col(position = position_dodge2()) +
      theme_prism(base_size = 10, axis_text_angle = 45) +
      scale_fill_manual(values = mycolorsLC[3:4]) +
      scale_y_continuous(guide = "prism_offset", expand = expansion(mult = c(0,NA))) +
      labs(y="Number of peaks", x="", title = "HBE day 7") +
      #geom_text(aes(label=n), position = position_dodge2(width = 0.9), vjust=-0.5, size=3) +
      theme(legend.position = c(0.75,0.9), aspect.ratio = 0.75) -> pNpeaks_HBEd7
    
    
    peaks_samples %>% filter(name %in% N27_pass) %>%
      filter(!(name %in% c(fiveUTRStartPeaks_N27mock,fiveUTRStartPeaks_N27cov2))) %>% 
      mutate(cond=if_else(grepl("cov2",cond),"SARS-CoV-2","Mock")) %>% 
      countPeaksPerBin(.) %>% ### Calculate counts per bins and normalize to spike-in
      mutate(countsPerBin=if_else(grepl("Mock",cond), countsPerBin * sf_N27mock, countsPerBin * sf_N27cov2 )) %>%
      mutate(type=if_else(dplyr::between(bin,0,2),"upstream",
                          if_else(dplyr::between(bin,2,3),"5'UTR",
                                  if_else(dplyr::between(bin,3,5.8),"CDS",
                                          if_else(dplyr::between(bin,5.9,6),"stop codon",
                                                  if_else(dplyr::between(bin,6.1,7),"3'UTR",
                                                          if_else(dplyr::between(bin,6.2,8),"downstream","NA"))))))) %>% 
      group_by(cond,type) %>% 
      mutate(n=round(sum(countsPerBin))) %>% 
      ungroup() %>%
      select(cond,type,n) %>% 
      distinct() %>% 
      filter(! grepl("stream",type) ) %>% 
      ggplot(aes(fct_relevel(type,order.type),n,fill=cond)) + 
      geom_col(position = position_dodge2()) +
      theme_prism(base_size = 10, axis_text_angle = 45) +
      scale_fill_manual(values = colors3[c(1,3)]) +
      scale_y_continuous(guide = "prism_offset", expand = expansion(mult = c(0,NA))) +
      labs(y="Number of peaks", x="", title = "HNE") +
      #geom_text(aes(label=n), position = position_dodge2(width = 0.9), vjust=-0.5, size=3) +
      theme(legend.position = c(0.75,0.9), aspect.ratio = 0.75) -> pNpeaks_N27
    
    # Vero
    countsPerBin_peaks_vero %>%
      mutate(type=if_else(dplyr::between(bin,0,2),"upstream",
                                                if_else(dplyr::between(bin,2,3),"5'UTR",
                                                        if_else(dplyr::between(bin,3,5.8),"CDS",
                                                                if_else(dplyr::between(bin,5.9,6),"stop codon",
                                                                        if_else(dplyr::between(bin,6.1,7),"3'UTR",
                                                                                if_else(dplyr::between(bin,6.2,8),"downstream","NA"))))))) %>% 
      group_by(cond,type) %>% 
      mutate(n=round(sum(countsPerBin))) %>% 
      ungroup() %>%
      select(cond,type,n) %>% 
      distinct() %>% 
      filter(!grepl("stream",type) ) %>% 
      ggplot(aes(fct_relevel(type,order.type),n,fill=cond)) + 
      geom_col(position = position_dodge2()) +
      theme_prism(base_size = 10, axis_text_angle = 45) +
      scale_fill_manual(values = mycolors) +
      scale_y_continuous(guide = "prism_offset", expand = expansion(mult = c(0,NA))) +
      labs(y="Number of peaks", x="", title = "Vero") +
      #geom_text(aes(label=n), position = position_dodge2(width = 0.9), vjust=-0.5, size=3) +
      theme(legend.position = c(0.75,0.9), aspect.ratio = 0.75) -> pNpeaks_vero
    
    
    pNpeaks <- (pNpeaks_HBEd4 | pNpeaks_HBEd7) / (pNpeaks_N27 | pNpeaks_vero)
    
    ggsave("/plots/numberOfPeaks_metagene_HBE_HNE_Vero.pdf", plot = pNpeaks, width = 10, height = 12, dpi = 300)
    
    #peaksPerBinPlot_HBEsamples + scale_x_continuous(breaks = seq(1,8,by=0.5)) + annotate(geom = "rect",xmin = 5.7, xmax = 6.1,ymin = 0, ymax = 200, alpha=0.5, fill="blue") + scale_y_continuous(breaks = seq(0,200,by=10))
    
    
#_________________________________________________________________________#
#---------- Select peaks in 20bp from start:---------------------------
#_________________________________________________________________________#


nearTSS <- function(df.gr,scaledPeaks,metagene_20FromStart_idx, cond){
  
  names20 <- join_overlap_inner_directed(df.gr %>% 
                                           mutate(peakStart=start,start=peakStart+peak.point.source,end=peakStart+peak.point.source+1,cond=cond),
                                         makeGRangesFromDataFrame(metagene_20FromStart_idx,keep.extra.columns = T)) %>% 
    as.data.frame() %>% select(name) %>% distinct() %>% unlist() 
   
  res <- scaledPeaks %>% filter(name %in% names20 & type=="cds" & dplyr::between(scaled_summit_bin,3,3.5))
  
  return(res)
}

distFromStart <- 20

metagene_50FromStart_chlsab_idx <- metagenes_chlsab_long %>% mutate(idx=row_number()) %>% filter(type=="cds") %>% mutate(cds_start=start,cds_end=end,cds_strand=strand) %>% 
  mutate(start=if_else(cds_strand=="+",cds_start,cds_end - distFromStart),
         end=if_else(cds_strand=="+",cds_start + distFromStart,cds_end)) #%>% dplyr::select(idx) %>% distinct()

metagene_50FromStart_human_idx <- metagenes_hs_long %>% mutate(idx=row_number()) %>% filter(type=="cds") %>% mutate(cds_start=start,cds_end=end,cds_strand=strand) %>% 
  mutate(start=if_else(cds_strand=="+",cds_start,cds_end - distFromStart),
         end=if_else(cds_strand=="+",cds_start + distFromStart,cds_end))


############## Vero -------------------------------


scaledPeaks_FromStart_veroNonInf <- nearTSS(vero_peaks.gr,scaledPeaks_vero,metagene_50FromStart_chlsab_idx, cond="Non-infected")

scaledPeaks_FromStart_B.1 <- nearTSS(wu_peaks.gr,scaledPeaks_vero,metagene_50FromStart_chlsab_idx, cond="B.1")

scaledPeaks_FromStart_B.1.1.7 <- nearTSS(uk_peaks.gr,scaledPeaks_vero, metagene_50FromStart_chlsab_idx,cond="B.1.1.7")

scaledPeaks_FromStart_B.1.351 <- nearTSS(sa_peaks.gr,scaledPeaks_vero, metagene_50FromStart_chlsab_idx,cond="B.1.351")

scaledPeaks_FromStart_veroNonInf %>% mutate(origStart=start,start=(origStart+peak.point.source)-4,end=(origStart +peak.point.source)+4) %>% select(seqnames,start,end,name,score,strand) %>% fwrite(.,"/genomeResearch_cov2_v2/tables/peaksNearTSS_vero_nonInfected.bed", col.names = F, sep="\t")

scaledPeaks_FromStart_veroNonInf %>% mutate(start=(start+peak.point.source)-10,end=start+10) %>% select(seqnames,start,end,name,score,strand) %>% fwrite(.,"/genomeResearch_cov2_v2/tables/peaksNearTSS_vero_nonInfected.bed", col.names = F, sep="\t")

################ HBE samples________________________

scaledPeaks_FromStart_LCmockD4 <- nearTSS(LC_mock_d4,scaledPeaks_LC,metagene_50FromStart_human_idx, cond="HBE mock day 4")

scaledPeaks_FromStart_LCmockD7 <- nearTSS(LC_mock_d7,scaledPeaks_LC,metagene_50FromStart_human_idx, cond="HBE mock day 7")
  #overlapPeakAndGenomeCoords2(LC_mock_d7 %>% mutate(cond="HBE Mock d7") %>% filter(name %in% cdsStartPeaks_LCmockd7$name),makeGRangesFromDataFrame(metagenes_hs_long,keep.extra.columns = T), tot.bins = tot.bins) %>% filter(type=="cds",scaled_summit_bin < 3.5)

scaledPeaks_FromStart_LCcov2D4 <- nearTSS(LC_cov2_d4,scaledPeaks_LC,metagene_50FromStart_human_idx, cond="HBE cov2 day 4")

scaledPeaks_FromStart_LCcov2D7 <- nearTSS(LC_cov2_d7,scaledPeaks_LC,metagene_50FromStart_human_idx, cond="HBE cov2 day 7")

scaledPeaks_20bpFromStart_LCmockD7 %>% countPeaksPerBin(.) %>% plotPeaksPerBinDensity()
scaledPeaks_20bpFromStart_LCmockD4 %>% countPeaksPerBin(.) %>% plotPeaksPerBinDensity()




#_____________ Plot peaks Without the first 20bp from TSS:_______________________________#

all20FromStart_<- c(scaledPeaks_20bpFromStart_LCmockD4$name,scaledPeaks_20bpFromStart_LCmockD7$name,scaledPeaks_20bpFromStart_LCcov2D4$name,scaledPeaks_20bpFromStart_LCcov2D7$name)

scaledPeaks_%>% filter(!(name %in% all20FromStart_LC)) %>% countPeaksPerBin(.) %>% plotPeaksPerBinDensity()


cdsStartPeaks_LCmockd7 %>% mutate(start=start-5,end=end+5) %>% dplyr::select(seqnames,start,end,score,name,strand) %>% fwrite(.,"/peakCalling/peaksCDSstart_LCmockD7.bed", col.names = F, sep="\t")



scaledPeaks_LCmockd7 %>% filter(name %in% cdsStartPeaks_LCmockd7$name)

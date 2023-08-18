
library(AnnotationDbi)
library(GenomicFeatures)
library(org.Hs.eg.db)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)


txdb_hg38 <- makeTxDbFromUCSC(genome='hg38', tablename='refGene')

peaksXiang <- fread("public_datasets/peaks_U2OS_GSE92867_WT.minus.UV.longest_xiang2017.bed", header = F, sep = "\t", col.names = c("refseqGene","start","end","name","score","strand"))
  

  peaksXiang <- peaksXiang %>% rowwise() %>% mutate(chr=select(org.Hs.eg.db,keys = refseqGene, keytype = "REFSEQ", columns = "CHR"))

  #left_join(peaksXiang, refseq_hg38, by=c("refseqGene"="refseq"))
  
  peaksXiang <- left_join(peaksXiang,transcripts(txdb_hg38) %>% as.data.frame(), by=c("refseqGene"="tx_name"))
  
  
  peaksXiang %>% mutate(start=start.x+start.y, end=end.x+end.y) %>% filter(seqnames %in% paste0("chr",c(seq(1,22),"X","Y"))) %>% mutate(score=".") %>% dplyr::select(seqnames,start,end,refseqGene,score,strand.x) %>% fwrite("public_datasets/peaks_U2OS_GSE92867_WT.minus.UV.longest_xiang2017_genomeCoords.bed", col.names = FALSE, sep="\t")
  

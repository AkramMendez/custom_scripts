
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)
library(Guitar)

txdb.hg38<-makeTxDbFromUCSC(genome="hg38", tablename="refGene")

guitarTxdb<-makeGuitarTxdb(txdb.hg38, txPrimaryOnly = FALSE)

control.bed<-"peak_summits_control.bed"

mettl3.bed<-"peak_summits_kd.bed"

mettl3.control.bed<-list(control.bed,mettl3.bed)


GuitarPlot(txTxdb = txdb.hg38,stBedFiles = mettl3.control.bed,headOrtail = T,enableCI = F,mapFilterTranscript = T,pltTxType = c("mrna"),stGroupName = c("Control","KD"))
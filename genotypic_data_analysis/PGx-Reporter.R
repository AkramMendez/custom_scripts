#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)
#How to run: Rscript --vanilla sillyScript.R data.to.pass
if(length(args)==0){
  stop("Provide the following parameters: \n args[1]= Path to Final-reports files (previously shrinked to the variants of interest), example:'../final_reports_by_EC/' \n args[2]= Path to Variants table, example: '../../../Cathalogs/variante.xlsx'\n args[3]= output path
", call.=FALSE)
}
library(readxl)
#args[1]= Path to Final-reports files (previously shrinked to the variants of interest), example: "../final_reports_by_EC/"
#args[2]= Path to Variants table, example: "../../../Cathalogs/variante.xlsx"
#args[3]= output path
#setwd("/home/akram/CD46/cd46_work/PGx/Data/")
#path<-"/home/akram/CD46/cd46_work/PGx/Data/final-reports/Final_reports_selected/"
#manifiesto<-read.delim(file="", header = T,stringsAsFactors = F, as.is = T, fileEncoding = "UTF-8")
#pgx.table<-read.csv(file = "", header = T,stringsAsFactors = F, as.is = T, fileEncoding = "UTF-8")
#final.report<-read.delim(file = "./", header = T,stringsAsFactors = F, as.is = T)
#pgx.table<-read_xlsx("../../../Cathalogs/PGx.xlsx",col_names = T)
#variant.table<-read_xlsx("../../../Cathalogs/variante.xlsx",col_names = T)
reports.path<-args[1]
#pgx.path<-args[2]
vars.path<-args[2]
path<-args[3]
files<-list.files(reports.path)
#pgx.table<-read_xlsx(pgx.path,col_names = T)
variant.table<-read_xlsx(vars.path,col_names = T)


#final.report<-read.delim(file = "../final_reports/final_reports_by_EC/EC1028.txt", header = T,stringsAsFactors = F, as.is = T, sep="\t")
sapply(files, function(file){
  
  cat(">File:", file, "\n")
#final.report<-read.delim(file = paste0("../final_reports_by_EC/",file), header = T,stringsAsFactors = F, as.is = T, sep="\t")
  final.report<-read.delim(file = paste0(reports.path,file), header = F,stringsAsFactors = F, as.is = T, sep="\t")
#if(ncol(final.report==5)){

colnames(final.report)<-c("Name","SampleID", "All1_TOP","All2_TOP", "All1-Fwd", "All2-Fwd","GC")
#colnames(final.report)<-c("Name","SampleID", "All1_TOP","All2_TOP", "GC_Score")
allele.columns<-c(3,4)
#}else{
#  colnames(final.report)<-c("Name","SampleID", "All1_TOP","All2_TOP", "GC_Score", "All1_Fwd","All2_Fwd")
#  allele.columns<-c(6,7)
#}


sample.name<-final.report$SampleID[1]

matched.name<-match(variant.table$Illumina_name, final.report$Name)

final.report.selected<-final.report[matched.name,]

get.names<-match(final.report.selected$Name, variant.table$Illumina_name)

final.report.selected<-cbind(variant.table[get.names,c("PK_Variant", "dbSNP_ID")],final.report.selected)
rownames(final.report.selected)<-NULL
#file<-paste0(path,"Selected_Variants_PGx_",sample.name,".txt")

write.table(final.report.selected,paste0(path,"Final-report_",sample.name,".txt"), col.names = T, fileEncoding = "UTF-8", sep="\t", row.names = F, append = F, quote = F)

#Clear environment
rm(list=ls())
})

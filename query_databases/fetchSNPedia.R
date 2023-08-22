#!/usr/bin/env Rscript
# Script for fetching SNP description data from SNPedia webpages
# The script retrieves genotypic and other related information in JSON format and process it to return a tabular file
args = commandArgs(trailingOnly=TRUE)
#How to run: Rscript --vanilla makeFinalReport.R Final-reports-folder variants-table-xlsx output-path
if(length(args)==0){
  stop("    Usage: Rscript --vanilla fetchSNPedia.R variants-list_(single column format) output-path \n", call.=FALSE)
}


library(SNPediaR)
library(jsonlite)
library(curl)
library(rjson)
library(RJSONIO)

# Function declaration

getGenotypeTags<-function(rsid,genotype){
  genotype.tag<-paste0(rsid,genotype)
  genotype.page<-getPages(genotype.tag)
  names(genotype.page)<-gsub("Rs","rs",names(genotype.page))
  cat("    Fetching genotype: ",genotype.tag, "\n")
  labels<-c("allele1","allele2","magnitude","repute","summary")
  if(!is.null(genotype.page[[1]])){
  tags<-t(extractTags(genotype.page[[genotype.tag]], tags = labels))  
  tags<-c(rsid,genotype,tags)
  }else{
    tags<-c(rsid,genotype,rep(0,length(labels)))
  }
  names(tags)<-c("rsid","genotype",labels)
  return(tags)
}

#Loading variants list
rsids<-read.table(file = args[1], col.names = F,stringsAsFactors = F, as.is = T, row.names = NULL)
output.path <- args[2]

names(rsids)<-"dbSNP_ID"

for(i in 1:length(rsids[[1]])){
rsid<-rsids$dbSNP_ID[i]
res<-getPages(rsid)
if(!is.null(unlist(res))){
snp.page<-t(sapply (res, extractSnpTags))
genotypes<-snp.page[,c("geno1","geno2","geno3")]
}else{genotypes<-NULL}

if(!is.null(genotypes)){
sapply(genotypes, function(i){
 cat("-",i,"\n")
 genotype.info<-t(getGenotypeTags(rsid,i))
  cat(genotype.info,"\n")
 write.table(genotype.info,paste0(output.path,"genotypes_from_SNPedia_variants_Manifest.v1.txt"),sep = "\t", quote = F, col.names = F, row.names = F, append = T)
 return(NULL)
})
}
}

#! /usr/bin/env Rscript
#   Developer: Akram Méndez
#   Script: updateSampleFiles.R
#   Description: Este script actualiza los valores de Sample_ID para nuevos Códigos de Barras, 
#   busca si existen en la base de datos y genera un archivo Samples y updateSamples.
#   Date: 19 Jun 2019
#   Código 46.

library(odbc)
library(DBI)
library(tidyverse)
library(data.table)
library(getPass)


args = commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  stop("    Usage: Rscript --vanilla updateSampleFiles.R Extraction_file.csv output_dir \n", call.=FALSE)
}

input<-args[1]
output<-path.expand(args[2])

if(is.na(output)){
  message("Setting output to current dir.")
  output<-"./"
}  

extraction<-fread(input,header = T,sep=";")


# Get Extraction file name
file_name<-unlist(input %>% str_split(.,"\\/")) %>% tail(1) %>% gsub("(\\.txt|\\.csv|\\.tsv)","",.,perl = T)

#Connect to C46_SamplesData database
#NOTE: The 'Driver' parameter must be changed depending on the installed driver, in AWS the driver is ODBC Driver 17 for SQL Server.
samplesData <- DBI::dbConnect(odbc::odbc(), Driver = "ODBC Driver 17 for SQL Server", Database = "C46_SamplesData", UID = getPass::getPass(msg = "Plase insert 'C46_SamplesData' \n USER NAME:",forcemask = T), PWD =getPass::getPass(msg = "C46_SamplesData \n PASSWORD:",forcemask = T), server = "codigo46.me", port = 1433)
#Assign default headers
extraction<-extraction[,c(2,3)]
colnames(extraction)<-c("Lab_label","Codigo_Barras")
extraction$Codigo_Barras<-as.character(extraction$Codigo_Barras)
barcodes<-extraction$Codigo_Barras
#Concatenate barcodes and query
string<-paste0(gsub("\\b","'",barcodes,perl = T),collapse = ",")  

#Search existing samples in the database
in_database<-dbGetQuery(samplesData,paste0("SELECT S.BarCode FROM Samples S WHERE S.BarCode IN (",string,") AND S.Sample_ID != ''" ))
samples<-semi_join(extraction,in_database,by=c("Codigo_Barras"="BarCode")) %>% select(Lab_label,Codigo_Barras)
colnames(samples)<-c("Sample_ID","BarCode")

#Get new samples in the extraction file
newSamples<-anti_join(extraction,in_database,by=c("Codigo_Barras"="BarCode")) %>% select(Lab_label,Codigo_Barras)

#Generate queries for updating samples
updateSamples<-sapply(1:nrow(newSamples), function(x){
  paste0("UPDATE Samples SET Sample_ID = '",newSamples[x,"Lab_label"],"' WHERE BarCode = '",newSamples[x,"Codigo_Barras"],"' AND Sample_ID = ''","\n","UPDATE EstadoSamples SET FK01_Estado = 6 WHERE BarCode = '",newSamples[x,"Codigo_Barras"],"'")
})


#Save Samples and updateSamples files:
message("Saving Samples and updateSamples files:")

if(nrow(samples)>0){
	message("Saving reprocessed samples for updating.")
	write_csv(samples,paste0(output,"Samples.",Sys.Date(),".csv"),col_names = T)
}else{
	message("There are no previously registered samples, the 'Samples' file have been not generated.")}
write.table(as.data.frame(updateSamples),paste0(output,"updateSamples_",file_name,"_",Sys.Date() %>% gsub("-","",.),".txt"),col.names = FALSE,row.names=FALSE,sep = "\t",quote =FALSE)

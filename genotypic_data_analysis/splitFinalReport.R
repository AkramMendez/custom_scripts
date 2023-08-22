#!/usr/bin/env Rscript
#Developer: Akram Méndez
#Script: splitFinalReport.R
#Description: Programa para separar un lote de FinalReports proveniente de GenomeStudio en muestras individuales.
#Copyright: Código 46, 2019.
#

args = commandArgs(trailingOnly=TRUE)

if(length(args)==0){
  stop("splitFinalReport.R: Este programa separa los datos crudos a partir de un archivo FinalReport generado en GenomeStudio para un lote de muestras.
  Por cada muestra genera un reporte de tipo 'SampleID_all_variants.txt' y los almacena en un directorio correspondiente a cada proyecto.
  Uso: Rscript --vanilla splitFinalReport.R FinalReport output_path(optional)
  <FinalReport>   Archivo FinalReport exportado desde GenomeStudio conteniendo los datos crudos para un lote de muestras.
  <output_path>   Ruta donde guardar los archivos separados *all_variants.txt para cada muestra.
                  Ruta por defaul: /media/datashare/data/Produccion/reports-all-variants/", call.=FALSE)
}

library(data.table)
library(tidyverse)
library(utils)

report_path <- path.expand(args[1])
output_path <- args[2]
#Extract Project name from FinalReport name (Anything before the "_FinalReport" string).
project_name<- str_split(report_path,pattern = "\\/", simplify=T) %>% last()
project_name<-str_split(project_name,"_FinalReport*", simplify = T) %>% first()

#check if output_path folder exists:


if(is.na(output_path)){
  output_path<-path.expand(paste0("./",project_name))
  if(dir.exists(output_path)){
  message("Setting output_path to (default):", output_path)
      }else{
          message("Making Project directory to default path:", output_path)
          dir.create(output_path)
          }
}else{
  output_path<-path.expand(paste(output_path,project_name,sep="/"))
    if(dir.exists(output_path)){
      message("Setting output_path to:", output_path)
    }else{
      message("Making Project directory to default path:", output_path)
      dir.create(output_path)
    }
  }

checkHeader<-function(header){
  #If header has 9 columns, check if it corresponds to an Fzero's FinalReport,
  #otherwise, check if it's a default FinalReport with 7 columns.
  	if(length(header)==5){
	message("        Processing FinalReport in Fzero's format (5 columns)")
    default_names<-c("SNP Name","Sample ID","Allele1 - Top","Allele2 - Top","GC Score")
}

if(length(header)==9){
    message("        Processing FinalReport in Fzero's format (9 columns)")
    default_names<-c("SNP Name","Sample ID","Allele1 - Top","Allele2 - Top","Allele1 - Forward","Allele2 - Forward","GC Score","B Allele Freq","Log R Ratio")
  }

if(length(header)==7){
    message("        Processing FinalReport with default headers (7 columns)")
  default_names<-c("SNP Name","Sample ID","Allele1 - Top","Allele2 - Top","Allele1 - Forward","Allele2 - Forward","GC Score")
}
  if(!all(header==default_names)){
    stop( "         Please provide a FinalReport with the 1st line with the correct header.\n",call.=FALSE)
    return(FALSE)
  }else{
    return(TRUE)
  }
}
#Read FinalReport
final_report<-fread(report_path, header = T, stringsAsFactors = F, sep="\t", fill = T, skip = 9, showProgress=T, colClasses=list(character=1:4))
#Get header:
header<-colnames(final_report)

#Check header, if correct, split FinalReport:
if(checkHeader(header)){
message("Reading FinalReport ...")
#Get Sample_IDs from FinalReport:
sample_ids<-unique(final_report$`Sample ID`)
# Split FinalReport:
message(paste0("Splitting FinalReports and saving to: ",output_path))
invisible(sapply(sample_ids,function(sample){
  fwrite(final_report[`Sample ID`==as.character(sample),],paste0(output_path,"/",as.character(sample),"_all_variants.txt"),sep = "\t")
  }))
message("Done.")
}else(stop("Incorrect header. Please provide a FinalReport with correct header.\n",call.=FALSE))

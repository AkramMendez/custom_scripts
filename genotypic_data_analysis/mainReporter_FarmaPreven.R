#!/usr/bin/env Rscript
#
# Developer:    Akram Méndez
# Script:       Process the raw final report batch to extract a file for each sample
#               containing all variants of the chip.
#               Código 46, 2018.
#

# Load data.table for fast processing of FinalReport tables
library(data.table)
library(utils)

# data.table requires also the library 'bit.64', install it with: install.packages('bit64')
args <- commandArgs(trailingOnly = TRUE)

if (length(args) == 0 ) {
  stop("         Farma-Preven reporter:		
         Main reporter to generate Farma and Prevención table reports. 
         Usage:                       Rscript --vanilla mainReporter.R reports-all-variants-dir output-dir(optional) \n
         <reports-all-variants-dir>  : Path to the folder containing individual reports or path to a single sample FinalReport.\n
         <output-dir>                : Path to save the FinalReport tables for Farma and Preven reports \n 
         (default: /media/datashare/data/EXPORT-TO-SERVER/FinalReports). \n
       ", call. = FALSE)
}


reports_all_variants <- args[1]


# Relative Paths:
# 1.  Path to the Variants table for Farma and Prevención reports.
# 2.  Path to the single sample FinalReports (previously extracted from batches with 'getSampleData.sh'.
# 3.  Path to the reports dir
#4.   Path to the sample_codes dir

# Input directory:
#variants_table_file <- "/media/datashare/data/Produccion/tables/v2/VariantsTable_produccion.txt"
variants_table_file <- "/home/akram/Documents/work/tablasProduccion/Variants.txt"
# Set Output directory:
if(is.na(args[2])){
  #output_dir<- "/media/datashare/data/Produccion/EXPORT-TO-SERVER/FinalReports"
  output_dir<- "/home/akram/Documents/work/"
  message("         Setting Output directory to: \n", output_dir)
}else{
  output_dir<- args[2]
  message("         Setting Output directory to: \n", output_dir)
}

checkHeader<-function(header){
if(length(header)==9){  
default_names<-c("SNP Name","Sample ID","Allele1 - Top","Allele2 - Top","Allele1 - Forward","Allele2 - Forward","GC Score","B Allele Freq","Log R Ratio")
}
if(length(header)==7){
default_names<-c("SNP Name","Sample ID","Allele1 - Top","Allele2 - Top","Allele1 - Forward","Allele2 - Forward","GC Score")
}
if(length(header)==5){
default_names<-c("SNP Name","Sample ID","Allele1 - Top","Allele2 - Top","GC Score")
}
  
if(!all(header==default_names)){
    stop( "         Please provide a FinalReport with the 1st line in the following column order:\n
         SNP Name, Sample ID, Allele1 - Top,Allele2 - Top, Allele1 - Forward, Allele2 - Forward, GC Score, B Allele Freq, Log R Ratio \n ")
    return(FALSE)
  }else{
    return(TRUE)
  }
}

saveToLogfile<-function(sample,project_name, ...){
  #logfile<-"/media/datashare/data/Produccion/logfiles/logfile_FinalReports_Farma_Preven.txt"
  date<-as.character(Sys.time())
if(file.exists(logfile)){
  message("         Saving sample ",sample, " in logfile.")
  write.table(t(c(sample,project_name,date)),file = logfile, col.names = F, fileEncoding = "UTF-8", sep = "\t", row.names = F, append = T, quote = F)
}else{  
  message("         Creating FinalReports logfile.")
  write.table(t(c(Sample_ID=sample,Batch=project_name,Date=date)),file = logfile, col.names = T, fileEncoding = "UTF-8", sep = "\t", row.names = F, append = F, quote = F)
}
}

matchSampleVariants <- function(report_path, project_name ,...) {
  # Load single sample data extracted from the raw FinalReports.
  message("         Loading file : \n", report_path)
  final_report <- fread(input = report_path, header = T, stringsAsFactors = F, sep = "\t", encoding = "UTF-8", blank.lines.skip = T,colClasses=list(character=1:4))
  if(checkHeader(colnames(final_report))){
    message("         Checking header ...")
	if(ncol(final_report)==9){
    	colnames(final_report) <- c("Name", "SampleID", "All1_TOP", "All2_TOP", "All1-Fwd", "All2-Fwd", "GC","B Allele Freq","Log R Ratio")
	  }
	if(ncol(final_report)==7){
    	colnames(final_report) <- c("Name", "SampleID", "All1_TOP", "All2_TOP", "All1-Fwd", "All2-Fwd", "GC")
  	}
	if(ncol(final_report)==5){
    	colnames(final_report) <- c("Name", "SampleID", "All1_TOP", "All2_TOP", "GC")
  	}
}
  # Get the sample name from the first field in the SampleID column
  if (length(unique(final_report$SampleID)) != 1) {
    stop("         Error: The file contains more than one SampleID.")
  }
  # Get the name for each sample as it's been processed
  sample_name <- final_report$SampleID[1]
  message("         Processing sample: ", sample_name, "\n")
  # Select FinalReport colums until the All1TOP and All2TOP:
  column_subset <- c("Name", "SampleID", "All1_TOP", "All2_TOP")
  final_report_subset <- final_report[, ..column_subset]
  # Match the Illumina name of variants in the sample report with Illumina Names in the CodigoFarma table Variants table.
  setkey(final_report_subset, Name)
  final_report_selected <- merge(variants_table, final_report_subset, by = "Name")
  # Reorder columns in the form of the FinalReports data structure:
  setcolorder(final_report_selected, c("FK01_Variant", "dbSNP_ID", "Name", "SampleID", "All1_TOP", "All2_TOP"))
  # Save all sample reports into a big single file (one file per project):
  # Name of the output file containing the data for the Farma and Preven report tables.
  output_file <- paste0(output_dir,"/FinalReports_",project_name,"_",date,".txt")
  # Check if output file exists, if no, create a new output file, if yes, overwrite the results
  # The output file are saved inside the output project dir with prefix "FinalReports_" and the project name.
  if (!file.exists(output_file)) {
    message("										Saving sample report:", sample_name, " \n 										in ", output_file)
    write.table(final_report_selected, output_file, col.names = T, fileEncoding = "UTF-8", sep = "\t", row.names = F, append = F, quote = F)
  }else {
    if (file.exists(output_file) && file.info(output_file)$size != 0) {
      message("         Saving sample ", sample_name)
      message("         Appending data to report: \n", output_file)
      write.table(final_report_selected, output_file, col.names = F, fileEncoding = "UTF-8", sep = "\t", row.names = F, append = T, quote = F)
    }
  }
  
  saveToLogfile(sample_name, project_name)
  return(as.character(sample_name))
  # Clear environment
  rm(list = ls())
}

# _______________Main Script_______________

#     0. Load and order the Variants table:
variants_table <- fread(input = variants_table_file, header = T, stringsAsFactors = F, sep = "\t", encoding = "UTF-8", blank.lines.skip = T)
colnames(variants_table)<-c("FK01_Variant","dbSNP_ID","Name")
setkey(variants_table, Name)

# Extract the project's name from the input argument, it can be a directory name (e.g. /Infinium_20180101/)
# or a single FinalReport name (e.g. FinalReport_Infinium_201808081.txt, EC1234_all_variants.txt)
project_name<-unlist(strsplit(reports_all_variants,"/"))
project_name<-project_name[length(project_name)]
project_name<-gsub("(FinalReport_|_all_variants.txt|.txt)","",project_name)
message("         Project: ",project_name)

date<-gsub("-","",Sys.Date())
samples_list<-NULL
# Read the single sample final reports:
if (file_test("-f", reports_all_variants)) {
  # This is a single final report file:
  message("         Processing single Final-report file: ", reports_all_variants, "\n")
  report_path<-normalizePath(reports_all_variants)
  sample_id<-matchSampleVariants(report_path, project_name, date)
  samples_list<-c(samples_list,sample_id)
}else {
  if (file_test("-d", reports_all_variants)) {
    # Final reports in the directory:
    message("         Processing all final reports in: ", reports_all_variants)
    files <- list.files(reports_all_variants)
    # Retrieve corresponding variants for each sample with "matchSampleVariants()" function.
    samples_list<-sapply(files, function(file){
      message("         Processing single sample FinalReport:", file, "\n")
      report_path<-normalizePath(paste0(reports_all_variants,file))
      sample_id<-matchSampleVariants(report_path, project_name, date)
      samples_list<-c(samples_list,sample_id)
      return(samples_list)
    })
  }
}



message("         Done.")

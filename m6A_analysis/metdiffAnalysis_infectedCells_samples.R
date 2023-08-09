#!/usr/bin/env Rscript

library(GenomicRanges)
library(GenomicAlignments)
library(rtracklayer)
library(IRanges)
library(tidyverse)
library(MeTDiff, lib.loc="/domus/h1/amendez/R/x86_64-pc-linux-gnu-library/4.0")


args = commandArgs(trailingOnly=TRUE)


gtf<-"/genomes/custom_reference_genomes/chlSab2_wuhCor1_ecoliK12_ncbiGenes_chroms_concat.gtf"
rev="/mappings/covid19_infected_cells/nodup_uniq_alns/rev"
fwd="/mappings/covid19_infected_cells/nodup_uniq_alns/fwd"
#inputdir<-"/mappings/covid19_viral_m6a_input/nodup_uniq_alns"
outdir<-"/peakCalling/metdiff/metdiff_infected_noninfected_cells_cov2"

#Input files Fwd:
input_WU_fwd=paste0(fwd,"/","COVID_19_Infected_cell_Input_1_S7_R1_001.chroms.nodup.uniq.sorted.fwd.bam")
input_SA_fwd=paste0(fwd,"/","COVID_SA_Infected_cell_Input_1_S9_R1_001.chroms.nodup.uniq.sorted.fwd.bam")
input_UK_fwd=paste0(fwd,"/","COVID_UK_Infected_cell_Input_1_S8_R1_001.chroms.nodup.uniq.sorted.fwd.bam")
input_vero_fwd=paste0(fwd,"/","Vero_Input_1_S13_R1_001.chroms.nodup.uniq.sorted.fwd.bam")

#IP files Fwd:
IP_WU_fwd=paste0(fwd,"/","COVID_19_Infected_cell_m6A_S10_R1_001.chroms.nodup.uniq.sorted.fwd.bam")
IP_SA_fwd=paste0(fwd,"/","COVID_SA_Infected_cell_m6A_S12_R1_001.chroms.nodup.uniq.sorted.fwd.bam")
IP_UK_fwd=paste0(fwd,"/","COVID_UK_Infected_cell_m6A_S11_R1_001_renamed.chroms.nodup.uniq.sorted.fwd.bam")
IP_vero_fwd=paste0(fwd,"/","Vero_m6A_S14_R1_001.chroms.nodup.uniq.sorted.fwd.bam")

#Names Fwd output file sufix:
name_WU_fwd="WU_infected_fwd"
name_SA_fwd="SA_infected_fwd"
name_UK_fwd="UK_infected_fwd"
names_vero_fwd="Vero_noninfected_fwd"

#Input files Rev:
input_WU_rev=paste0(rev,"/","COVID_19_Infected_cell_Input_1_S7_R1_001.chroms.nodup.uniq.sorted.rev.bam")
input_SA_rev=paste0(rev,"/","COVID_SA_Infected_cell_Input_1_S9_R1_001.chroms.nodup.uniq.sorted.rev.bam")
input_UK_rev=paste0(rev,"/","COVID_UK_Infected_cell_Input_1_S8_R1_001.chroms.nodup.uniq.sorted.rev.bam")
input_vero_rev=paste0(rev,"/","Vero_Input_1_S13_R1_001.chroms.nodup.uniq.sorted.rev.bam")

IP_WU_rev=paste0(rev,"/","COVID_19_Infected_cell_m6A_S10_R1_001.chroms.nodup.uniq.sorted.rev.bam")
IP_SA_rev=paste0(rev,"/","COVID_SA_Infected_cell_m6A_S12_R1_001.chroms.nodup.uniq.sorted.rev.bam")
IP_UK_rev=paste0(rev,"/","COVID_UK_Infected_cell_m6A_S11_R1_001_renamed.chroms.nodup.uniq.sorted.rev.bam")
IP_vero_rev=paste0(rev,"/","Vero_m6A_S14_R1_001.chroms.nodup.uniq.sorted.rev.bam")

name_WU_rev="WU_infected_rev"
name_SA_rev="SA_infected_rev"
name_UK_rev="UK_infected_rev"
names_vero_rev="Vero_noninfected_rev"

setwd(outdir)



IP_samples_rev<-c(IP_WU_rev=IP_WU_rev ,IP_SA_rev= IP_SA_rev,IP_UK_rev=IP_UK_rev,IP_vero_rev=IP_vero_rev)
input_samples_rev<-c(input_WU_rev=input_WU_rev,input_SA_rev=input_SA_rev,input_UK_rev=input_UK_rev,input_vero_rev=input_vero_rev)

IP_samples_fwd<-c(IP_WU_fwd=IP_WU_fwd,IP_SA_fwd=IP_SA_fwd,IP_UK_fwd=IP_UK_fwd,IP_vero_fwd=IP_vero_fwd)
input_samples_fwd<-c(input_WU_fwd=input_WU_fwd,input_SA_fwd=input_SA_fwd,input_UK_fwd=input_UK_fwd,input_vero_fwd=input_vero_fwd)

v1<-v2<-c("WU","SA","UK","Vero")

IP_file_comparisions_rev<-expand.grid(paste0("IP_",v1,"_rev"),paste0("IP_",v2,"_rev"))
input_file_comparisions_rev<-expand.grid(paste0("IP_",v1,"_rev"),paste0("IP_",v2,"_rev"))

colnames(IP_file_comparisions_rev)<-c("reference","contrast")
colnames(input_file_comparisions_rev)<-c("reference","contrast")

IP_file_comparisions_rev %>% filter(reference!=contrast) -> IP_file_comparisions_rev
input_file_comparisions_rev %>% filter(reference!=contrast) -> input_file_comparisions_rev


IP_file_comparisions_fwd<-expand.grid(paste0("IP_",v1,"_fwd"),paste0("IP_",v2,"_fwd"))
input_file_comparisions_fwd<-expand.grid(paste0("IP_",v1,"_fwd"),paste0("IP_",v2,"_fwd"))

colnames(IP_file_comparisions_fwd)<-c("reference","contrast")
colnames(input_file_comparisions_fwd)<-c("reference","contrast")

IP_file_comparisions_fwd %>% filter(reference!=contrast) -> IP_file_comparisions_fwd
input_file_comparisions_fwd %>% filter(reference!=contrast) -> input_file_comparisions_fwd

message("Starting strain comparisons for Rev strand data:")

for(i in 1:nrow(IP_file_comparisions_rev)){ 
    IP_ref<-IP_file_comparisions_rev[i,"reference"]
    input_ref<-input_file_comparisions_rev[i,"reference"]
    IP_contrast<-IP_file_comparisions_rev[i,"contrast"]
    input_contrast<-input_file_comparisions_rev[i,"contrast"]
   
comparison_name<-gsub("IP_","",paste(IP_file_comparisions_rev[i,"reference"],"vs",IP_file_comparisions_rev[i,"contrast"],sep = "_"))   

cat("Comparing strain samples ",IP_samples_rev[IP_ref]," and ",input_samples_rev[input_ref],"\n")

cat("Against ",IP_samples_rev[IP_contrast]," and ",input_samples_rev[input_contrast],"\n")

    metdiff(GENE_ANNO_GTF = gtf ,IP_BAM = IP_samples_rev[IP_ref],INPUT_BAM = input_samples_rev[input_ref], TREATED_IP_BAM = IP_samples_rev[IP_contrast] ,TREATED_INPUT_BAM = input_samples_rev[input_contrast], EXPERIMENT_NAME = comparison_name)

}


message("Starting strain comparisons for Fwd strand data:")
for(i in 1:nrow(IP_file_comparisions_fwd)){ 
    
    IP_ref<-IP_file_comparisions_fwd[i,"reference"]
    input_ref<-input_file_comparisions_fwd[i,"reference"]
    IP_contrast<-IP_file_comparisions_fwd[i,"contrast"]
    input_contrast<-input_file_comparisions_fwd[i,"contrast"]
   
    comparison_name<-gsub("IP_","",paste(IP_file_comparisions_fwd[i,"reference"],"vs",IP_file_comparisions_fwd[i,"contrast"],sep = "_"))   

    cat("Comparing strain samples ",IP_samples_fwd[IP_ref]," and ",input_samples_fwd[input_ref],"\n")

    cat("Against ",IP_samples_fwd[IP_contrast]," and ",input_samples_fwd[input_contrast],"\n")

    metdiff(GENE_ANNO_GTF = gtf ,IP_BAM = IP_samples_fwd[IP_ref],INPUT_BAM = input_samples_fwd[input_ref], TREATED_IP_BAM = IP_samples_fwd[IP_contrast] ,TREATED_INPUT_BAM = input_samples_fwd[input_contrast], EXPERIMENT_NAME = comparison_name)

}

message("Done.")
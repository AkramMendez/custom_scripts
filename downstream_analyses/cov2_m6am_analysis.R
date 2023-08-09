library(plyranges)
library(rtracklayer)

list.files("/peakCalling/peakCallingViralGenome/peaksInfectedCells_mappedTo_ViralGenome", pattern = "*.narrowPeak", recursive = T, full.names = T)


cov2_uk_fwd <- rtracklayer::import("/peakCalling/peakCallingViralGenome/peaksInfectedCells_mappedTo_ViralGenome/peaks_nodup_uniq_UK_infectedMappedToCov2_fwd/peaks_nodup_uniq_UK_infectedMappedToCov2_fwd_peaks.narrowPeak")

cov2_uk_rev <- rtracklayer::import("/peakCalling/peakCallingViralGenome/peaksInfectedCells_mappedTo_ViralGenome/peaks_nodup_uniq_UK_infectedMappedToCov2_rev/peaks_nodup_uniq_UK_infectedMappedToCov2_rev_peaks.narrowPeak")


cov2_sa_fwd <- rtracklayer::import("/peakCalling/peakCallingViralGenome/peaksInfectedCells_mappedTo_ViralGenome/peaks_nodup_uniq_SA_infectedMappedToCov2_fwd/peaks_nodup_uniq_SA_infectedMappedToCov2_fwd_peaks.narrowPeak")

cov2_sa_rev <- rtracklayer::import("/peakCalling/peakCallingViralGenome/peaksInfectedCells_mappedTo_ViralGenome/peaks_nodup_uniq_SA_infectedMappedToCov2_rev/peaks_nodup_uniq_SA_infectedMappedToCov2_rev_peaks.narrowPeak")

cov2_wu_fwd <- rtracklayer::import("/peakCalling/peakCallingViralGenome/peaksInfectedCells_mappedTo_ViralGenome/peaks_nodup_uniq_WU_infectedMappedToCov2_fwd/peaks_nodup_uniq_WU_infectedMappedToCov2_fwd_peaks.narrowPeak")

cov2_wu_rev <- rtracklayer::import("/peakCalling/peakCallingViralGenome/peaksInfectedCells_mappedTo_ViralGenome/peaks_nodup_uniq_WU_infectedMappedToCov2_rev/peaks_nodup_uniq_WU_infectedMappedToCov2_rev_peaks.narrowPeak")

cov2_genes_gtf <- rtracklayer::import(gtf_cov2)

cov2_genes_gtf_20bp <- cov2_genes_gtf %>% as.data.frame() %>% filter(type=="transcript") %>% mutate(startOrig=start,endOrig=end) %>% mutate(start=start,end=start+20) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)

join_overlap_inner_directed(cov2_uk_fwd %>% as.data.frame() %>% mutate(strand=if_else(grepl("_fwd",name),"+","-")) %>% mutate(origStart=start,origEnd=end) %>% mutate(start=origStart+peak,end=origStart+peak) %>% makeGRangesFromDataFrame(keep.extra.columns = T), cov2_genes_gtf_20bp)

join_overlap_inner_directed(cov2_uk_rev %>% as.data.frame() %>% mutate(strand=if_else(grepl("_fwd",name),"+","-")) %>% mutate(origStart=start,origEnd=end) %>% mutate(start=origStart+peak,end=origStart+peak) %>% makeGRangesFromDataFrame(keep.extra.columns = T), cov2_genes_gtf_20bp)


join_overlap_inner_directed(cov2_sa_fwd %>% as.data.frame() %>% mutate(strand=if_else(grepl("_fwd",name),"+","-")) %>% mutate(origStart=start,origEnd=end) %>% mutate(start=origStart+peak,end=origStart+peak) %>% makeGRangesFromDataFrame(keep.extra.columns = T), cov2_genes_gtf_20bp)

join_overlap_inner_directed(cov2_sa_rev %>% as.data.frame() %>% mutate(strand=if_else(grepl("_fwd",name),"+","-")) %>% mutate(origStart=start,origEnd=end) %>% mutate(start=origStart+peak,end=origStart+peak) %>% makeGRangesFromDataFrame(keep.extra.columns = T), cov2_genes_gtf_20bp)


join_overlap_inner_directed(cov2_wu_fwd %>% as.data.frame() %>% mutate(strand=if_else(grepl("_fwd",name),"+","-")) %>% mutate(origStart=start,origEnd=end) %>% mutate(start=origStart+peak,end=origStart+peak) %>% makeGRangesFromDataFrame(keep.extra.columns = T), cov2_genes_gtf_20bp)

join_overlap_inner_directed(cov2_wu_rev %>% as.data.frame() %>% mutate(strand=if_else(grepl("_fwd",name),"+","-")) %>% mutate(origStart=start,origEnd=end) %>% mutate(start=origStart+peak,end=origStart+peak) %>% makeGRangesFromDataFrame(keep.extra.columns = T), cov2_genes_gtf_20bp)

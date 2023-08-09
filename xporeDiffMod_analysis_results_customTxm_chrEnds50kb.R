library(tidyverse)
library(data.table)
library(plyranges)
library(rtracklayer)


chrEndCoords <- rtracklayer::import("~/gitrepos/bioinfopipes/lab/CHM13_hg38chrY_chromEndCoordinates_50kb.genome",format="BED")

chrEndCoords <- chrEndCoords %>% as.data.frame() %>% mutate(strand=if_else(end > 50000,"+","-")) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T)


# Load xPore diffmod table resutls, for a description of each column check: https://xpore.readthedocs.io/en/latest/outputtable.html#outputtable
xptable_custom <- fread("~/MondalLab/terra/nanoporeTerra/xpore_diffmod/xpore_diffmod_customTxm/diffmod.table", sep = ",",header = T)


# Translate custom assembly coordinates to genomic coordinates:
xptable_custom <- xptable_custom %>%
  separate(id,
           into = c("arm", "range"),
           sep = "\\:\\:") %>%
  mutate(range = gsub("chr\\d+\\:", "", range)) %>%
  separate(range,
           into = c("chr_start", "chr_end"),
           sep = "-") %>%
  mutate(
    chr = gsub("\\_.*","", arm),
    start = (as.numeric(chr_start) + position) - 2,
    end = (as.numeric(chr_start) + position) + 3
  )

colnames(xptable_custom) <- gsub("-",".",colnames(xptable_custom))

DRACH <- "(A|G|T)(A|G)AC(A|C|T)"
NNACN <- "..AC."
NNANN <- "..A.."
nonCanon<- c("GAACU","UAAAG","GAAGA","GAAUU")

xptable_custom %>% filter(
  diff_mod_rate_CTRLSH_vs_METTL3KD > 0,
  grepl(NNANN, kmer, perl = T),
  pval_CTRLSH_vs_METTL3KD < 0.05
) %>% mutate(strand=".") %>% 
  dplyr::select(chr, start, end, kmer, position, strand) %>% 
  fwrite(
  .,
  "~/MondalLab/terra/nanoporeTerra/xpore_diffmod/xpore_diffmod_customTxm/coordinates_tx_to_genome_Custom_xpore_chrEnds50kb_NNANN.bed", col.names = F,sep = "\t"
)



# Compare peaks obtained from every assembly analysis woith xPore
modified_custom <- xptable_custom %>% filter(
  diff_mod_rate_CTRLSH_vs_METTL3KD > 0,
  grepl(NNANN, kmer, perl = T),
  pval_CTRLSH_vs_METTL3KD < 0.05 ,
  mod_rate_CTRLSH.REP1 > 0.25
) %>% mutate(strand = ".") %>% 
  dplyr::select(chr, start, end, kmer, position, strand,arm) %>% 
  mutate(start=start-150,end=end+150) %>% # Extend sites +-150bp around
  makeGRangesFromDataFrame(.,keep.extra.columns = T)

# xPore sites passing criteria from StringTie assembly on Nanopore Ctrl-sh data
modified_stringTieNano <- passing_sites_genomicCoors_str %>%
  as.data.frame() %>% mutate(strand == ".") %>%
  rename(start = coord_start, end = coord_end) %>%
  mutate(start = start - 150, end = end + 150) %>% # Extend sites +-150bp around
  makeGRangesFromDataFrame(., keep.extra.columns = T)

modified_stringtie_shortReads <-
  passing_sites_genomicCoors_str_shorReads %>%
  as.data.frame() %>% mutate(strand == ".") %>%
  rename(start = coord_start, end = coord_end) %>%
  mutate(start = start - 150, end = end + 150) %>% # Extend sites +-150bp around
  makeGRangesFromDataFrame(., keep.extra.columns = T)

join_overlap_inner(modified_custom,peaksCsh) # Overlapping sites custom assembly 50kb and Ctrl-sh m6A-RIP peaks

join_overlap_inner(modified_stringTieNano,peaksCsh) # Overlapping sites Stringtie Nanopore assembly 50kb and Ctrl-sh m6A-RIP peaks

join_overlap_inner(modified_stringtie_shortReads,peaksCsh) # Overlapping Stringtie Short-Reads and Ctrl-sh m6A-RIP peaks

join_overlap_inner(modified_stringtie_shortReads %>% as.data.frame() %>% mutate(id=row_number()) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),peaksCsh) %>% as.data.frame() %>% dplyr::select(id) %>% unlist() -> ids_shortReads

join_overlap_inner(modified_stringTieNano %>% as.data.frame() %>% mutate(id=row_number()) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),peaksCsh) %>% as.data.frame() %>% dplyr::select(id) %>% unlist() -> ids_longReads

modified_stringtie_shortReads %>% as.data.frame() %>% mutate(id=row_number()) %>% filter(id %in% ids_shortReads) %>% mutate(start=start+150,end=end-150,range=paste0(seqnames,":",start,"-",end))

modified_stringTieNano %>% as.data.frame() %>% mutate(id=row_number()) %>% filter(id %in% ids_longReads) %>% mutate(start=start+150,end=end-150,range=paste0(seqnames,":",start,"-",end))


modified_stringtie_shortReads %>% as.data.frame() %>% mutate(id=row_number()) %>% filter(id %in% ids_shortReads)

join_overlap_inner(
  join_overlap_inner(
    join_overlap_inner(modified_stringTieNano, peaksCsh),
    join_overlap_inner(modified_custom, peaksCsh)
  ),
  modified_stringtie_shortReads
)


join_overlap_inner(modified_stringtie_shortReads,peaksCsh) %>% filter(grepl("GAACT|TAAAG|GAAGA|GAATT|(A|G|T)(A|G)AC(A|C|T)", kmer))

join_overlap_inner(modified_custom,peaksCsh) %>% filter(grepl("GAACT|TAAAG|GAAGA|GAATT|(A|G|T)(A|G)AC(A|C|T)", kmer))

# Overlap m6a-RIP and Nanopores Txm
join_overlap_inner(passing_sites_genomicCoors_str %>% mutate(start=start-150,end=end+150) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),
                   mf1_tumor_peaks)

# Overlap m6a-RIP and Short-reads Txm
join_overlap_inner(passing_sites_genomicCoors_str_shorReads %>% mutate(start=start-150,end=end+150) %>% makeGRangesFromDataFrame(.,keep.extra.columns = T),
                   mf1_tumor_peaks
)



##########-------------------------------------------Sequence Logos----------------------------------------------------
#### Generate sequence logos from predicted modified k-mers: (https://omarwagih.github.io/ggseqlogo/)
library(ggseqlogo)
library(patchwork)
#Differentially modified k-mers, xPore on Stringtie's assembly using Nanopore reads:
kmers_str <- xptable_str %>% filter(diff_mod_rate_CTRLSH_vs_METTL3KD >  0, pval_CTRLSH_vs_METTL3KD < 0.05, grepl("..A..",kmer)) %>% dplyr::select(kmer) %>% unlist()
kmers_str_shortReads <- xptable_str_shorReads %>% filter(diff_mod_rate_CTRLSH_vs_METTL3KD >  0, pval_CTRLSH_vs_METTL3KD < 0.05, grepl("..A..",kmer)) %>% dplyr::select(kmer) %>% unlist()

kmers_str_custom <- xptable_custom %>% filter(diff_mod_rate_CTRLSH_vs_METTL3KD >  0, pval_CTRLSH_vs_METTL3KD < 0.05, grepl("..A..",kmer)) %>% dplyr::select(kmer) %>% unlist()

ggseqlogo(kmers_str, method="prob") + labs(title="Nanopore assembly") -> logoKmers_str
ggseqlogo(kmers_str_shortReads, method="prob") + labs(title="Short-reads assembly") -> logoKmers_str_shortReads
ggseqlogo(kmers_str_custom, method="prob") + labs(title="Custom assembly") -> logoKmers_str_custom

#ggseqlogo(c(kmers_str,kmers_str_shortReads), method="prob")

logoKmers_str + logoKmers_str_shortReads + logoKmers_str_custom -> p
#ggsave("~/MondalLab/terra/plots/sequenceLogo_kmers_diffMod0.25_pval0.05_xPore_custom_nanopore_shortReads_assemblies_NNANN_context.pdf", plot = p, width = 12, height = 4)
#ggsave("~/MondalLab/terra/plots/sequenceLogo_kmers_diffMod0.25_pval0.05_xPore_custom_nanopore_shortReads_assemblies.pdf", plot = p, width = 12, height = 4)
ggsave("~/MondalLab/terra/plots/sequenceLogo_kmers_diffModAbove0_pval0.05_xPore_custom_nanopore_shortReads_assemblies_NNANN_context.pdf", plot = p, width = 12, height = 4)
##########-----------------------------------------------------------------------------------------------



##### Generate sequence logos from predicted modified k-mers genome-wide

xptable_genwide_t2t <- fread("~/MondalLab/terra/nanoporeTerra/xpore_diffmod/diffmod.table", sep = ",",header = T)

xptable_genwide_t2t %>% filter( diff_mod_rate_CTRLSH_vs_METTL3KD >  0, 
                                pval_CTRLSH_vs_METTL3KD < 0.05,
                                `mod_rate_CTRLSH-REP1` > 0.9,
                                grepl("..A..",kmer)) %>% 
  dplyr::select(kmer) %>% group_by(kmer) %>% 
  dplyr::count() %>% 
  ungroup() %>% 
  mutate(total=sum(n)) %>% 
  mutate(perc=(n/total)*100) %>% 
  arrange(desc(perc)) %>%
  top_n(100,wt = perc) %>% 
  mutate(drach= if_else(grepl("(A|G|T)(A|G)AC(A|C|T)",kmer,perl = T),"DRACH","Non-DRACH")) %>% 
  ggplot(aes(reorder(kmer,desc(n)),n,fill=drach)) + geom_col() + scale_fill_manual(values = c("steelblue","grey70")) + 
  theme_prism(axis_text_angle = 45) + scale_y_continuous(expand = expansion(mult = c(0,0.01))) + 
  labs(x="k-mer",y="Frequency", title = "Top 100 k-mers \n xPore T2T genome-wide NNACN")

xptable_genwide_t2t %>% filter( diff_mod_rate_CTRLSH_vs_METTL3KD >  0 ,pval_CTRLSH_vs_METTL3KD < 0.05, grepl("..AC.",kmer)) %>% ggplot(aes(kmer,`mod_rate_CTRLSH-REP1`)) + geom_boxplot()


kmers_t2t <- xptable_genwide_t2t %>% filter( diff_mod_rate_CTRLSH_vs_METTL3KD >  0,
                                             pval_CTRLSH_vs_METTL3KD < 0.05, 
                                             grepl("..AC.",kmer)) %>% dplyr::select(kmer) %>% unlist()


ggseqlogo(kmers_t2t, method="prob") + labs(title="Nanopore T2T") #-> logoKmers_t2t

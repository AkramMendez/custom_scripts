library(tidyverse)

gseanames <- c("GENE_SET","RANK_IN_GENE_SET","SYMBOL","RANK_IN_GENE_LIST","RANK_METRIC_SCORE","RUNNING_ES","CORE_ENRICHMENT","")

gsea_kd_vs_ctrl_Untreated_cellType1 <- fread("/GSEA/kd_vs_ctrl_Untreated_cellType1_MigDBselected.GseaPreranked.1660825992554/gsea_tables/gsea_gathered_kd_vs_ctrl_Untreated_cellType1_MigDBselected.tsv", sep = "\t",header = F, col.names = gseanames)[,c(1:7)]

gsea_kd_vs_ctrl_treated_cellType1 <- fread("/GSEA/kd_vs_ctrl_treated_cellType1_MigDBselected.GseaPreranked.1660825426626/gsea_tables/gsea_gathered_kd_vs_ctrl_treated_cellType1_MigDBselected.tsv", sep = "\t",header = F, col.names = gseanames)[,c(1:7)]

gsea_kd_vs_ctrl__cellType2 <- fread("/GSEA/kd_vs_ctrl__cellType2_MigDBselected.GseaPreranked.1660825731902/gsea_tables/gsea_gathered_kd_vs_ctrl__cellType2_MigDBselected.tsv", sep = "\t",header = F, col.names = gseanames)[,c(1:7)]

left_join(
  gsea_kd_vs_ctrl_Untreated_cellType1 %>% mutate(significant_in_diffExpAnalysis=if_else(SYMBOL %in% unique(signif_cellType1_untreated_kdvsctrl$gene), "YES","NO")), signif_cellType1_untreated_kdvsctrl %>% dplyr::select(gene,log2FoldChange,padj), 
  by=c("SYMBOL"="gene")) %>% openxlsx::write.xlsx(.,"/GSEA/gatheredResults_GSEA_kd_vs_ctrl_Untreated_cellType1_MigDBselected.xlsx", overwrite = T)


left_join(
  gsea_kd_vs_ctrl_treated_cellType1 %>% mutate(significant_in_diffExpAnalysis=if_else(SYMBOL %in% unique(signif_cellType1_treated_kdvsctrl$gene), "YES","NO")), signif_cellType1_treated_kdvsctrl %>% dplyr::select(gene,log2FoldChange,padj), 
  by=c("SYMBOL"="gene")) %>% openxlsx::write.xlsx(.,"/GSEA/gatheredResults_GSEA_kd_vs_ctrl_treated_cellType1_MigDBselected.xlsx", overwrite = T)


left_join(
  gsea_kd_vs_ctrl__cellType2 %>% mutate(significant_in_diffExpAnalysis=if_else(SYMBOL %in% unique(signif__cellType2_kdvsctrl$gene), "YES","NO")), signif__cellType2_kdvsctrl %>% dplyr::select(gene,log2FoldChange,padj), 
  by=c("SYMBOL"="gene")) %>% openxlsx::write.xlsx(.,"/GSEA/gatheredResults_GSEA_kd_vs_ctrl_cellType2_MigDBselected.xlsx", overwrite = T)

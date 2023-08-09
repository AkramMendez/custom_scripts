library(tximport)
library(tximeta)
library(DESeq2)
library(tidyverse)
library(data.table)
library(openxlsx)

# Path to salmon quantification folders:
salmonDir <- "/salmon_quant/salmon_quant_cellType1_cellType2_batch/"

#Get complete paths for Salmon quantification files at the transcript level, Salmon files are named as "quant.sf"
salmon_samples <- list.files(salmonDir,recursive = T,pattern = "quant.sf",full.names = T)

# Read sample information table, it can be a CSV file: 
# This table must have at least the following columns so as to Tximport to read the files: "files", "names", "replicate", "condition"
metadata_salmon <- fread("salmonSamplesInfo_coldata_cellType1.csv", sep=",", header = T)

#Import salmon files and summarize counts at the gene level:
library(org.Sc.sgd.db) # Saccharomyces cerevicie annotations.You can install it with: BiocManager::install("org.Sc.sgd.db", update=FALSE)

    salmon_table <- tximeta(metadata_salmon,type = "salmon", skipSeqinfo = T) # tximeta will use the org.Sc.sgd.db loaded in the environment.
    salmon_table <- addIds(salmon_table_be2,"SYMBOL",gene = T)
    geneLevelSalmon <- summarizeToGene(salmon_table)


# Make a DESeqDataSet object from imported Salmon counts:
    ddsSalmon <- DESeqDataSet(geneLevelSalmon, design = ~ condition) #the 'design' parameter takes a formula to make several comparisons, in this case we wan't to compare two conditions, i.e. treated vs untreated, knockdown vs control

# Reassign 'factor' levels to reference conditions, this is important to set a correct order of the conditions we want to compare such that DESeq can understand which is the reference condition and which is the treatment condition:
    ddsSalmon$condition <- factor(ddsSalmon$condition, c("control","knockdown"))

# Run DESeq2 for differential expression analysis using the ddsSalmon object that contains metadata, counts and 'design' matrix for desired comparisons:
    ddsSalmon <- DESeq(ddsSalmon,modelMatrixType = "standard", fitType = "local")

# After running DESeq2, we can check the names of the results to make sure it used the correct comparisions, 
# for example it can give us something like this: "type_knowkcdow_vs_control", 
# meaning that it compared the counts in the knockdown (numerator) against  control counts (denominator).
# Therefore, genes with a positive log2FoldChange will be increased in knowckdown compared to the control condition, and negative log2FoldChanges will be genes lowly expressed after knockdown as compared to the control condition.

# Extract Differential expression results
    results_diffExpression <- results(ddsSalmon) %>% as.data.frame() %>% mutate(contrast="knockdown_vs_control") # We extract the log2 Foldchanges and add a column 'contrast' to keep track of the comparison we already made.

# We can apply an extra correction step that adjusts the variations in log2FoldChanges to reduce the number of false-positives that show high dispersion or low-expression:
# Note that the lfcShrink 'coef' parameter takes the contrasted condition name we got from the Salmon object, in this case it took a vector names as "condition_knowkdown_control"
# The type of fold change adjustmen used is "apeglm" which is a sensitive method for adjusting the fold changes:

    results_adjusted_LFC <- lfcShrink(ddsSalmon,coef = "condition_knockdown_vs_control",type = 'apeglm') %>% as.data.frame() %>% mutate(contrast="knockdown_vs_control") # We extract the adjusted log2 Foldchanges and add a column of the comparison we already made


#### Saving the results:

# Saving original differential expression analysis (adjusted by sample size, gene dispersion, and adjusted p-values):
openxlsx::write.xlsx(results_diffExpression,"/data/tables/differentialExpression_deseq2_knowkdown_vs_control.xls")

# Saving differential expression analysis with adjusted fold changes (discarded genes with low-expression and high variation):
openxlsx::write.xlsx(results_adjusted_LFC,"/data/tables/differentialExpression_deseq2_adjustedFoldChanges_knowkdown_vs_control.xls")


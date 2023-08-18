library(tidyverse)
library(GenomicFeatures)

# Set up filenames. Only change this line to point to the GTF.
gtf <- "gencode.v28.annotation.gtf"

# These filenames automatically generated from the above.
txdb.filename <- gsub("gtf", "sqlite", gtf)
tx2gene.filename <- gsub("gtf", "tx2gene.csv", gtf)


# Check to see if the txdb database file exists.
if (!file.exists(txdb.filename)) {
  ## If not, make it. Do it once. This takes a while.
  message(paste(txdb.filename, "doesn't exist. Creating and saving to disk..."))
  txdb <- makeTxDbFromGFF(gtf)
  saveDb(txdb, txdb.filename)
  message("Done.")
} else {
  ## If it already exists, load it from file, quickly!
  message(paste(txdb.filename, "found. Loading..."))
  txdb <- loadDb(txdb.filename)
  message("Done.")
}

keytypes(txdb)

# Create the tx2gene
tx2gene <- mapIds(txdb,
                  keys=keys(txdb, "GENEID"),
                  column="TXNAME",
                  keytype="GENEID",
                  multiVals="list") %>%
  enframe(name="ensgene", value="enstxp") %>%
  unnest() %>%
  dplyr::select(enstxp, ensgene) %>%
  distinct()

# Write out tx2gene with the header
write_csv(tx2gene, tx2gene.filename)

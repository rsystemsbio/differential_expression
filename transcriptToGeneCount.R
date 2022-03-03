#This script merges transcript counts (Kallisto Output) into gene counts
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("EnsDb.Hsapiens.v86")
BiocManager::install("tximportData")
BiocManager::install("tximport")
library(tximportData)
library(tximport)
kallisto_dir = ""
metadata_file = 'csample.txt'
metadata <- read.table(metadata_file, header = TRUE)
metadata

#BiocManager::install("TxDb.Hsapiens.UCSC.hg38.knownGene")
library(BSgenome.Hsapiens.UCSC.hg38)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
txdb <- BSgenome.Hsapiens.UCSC.hg38
columns(txdb)
txdb
keytypes(txdb)
k <- keys(txdb, keytype = "GENEID")
k
df <- select(txdb, keys=k, columns="TXNAME", keytype="GENEID")
df
#df <- select(txdb, keys = k, keytype = "GENEID", columns = "TXNAME")
tx2gene <- df[, 2:1]

files <- file.path(kallisto_dir, metadata$run, "abundance.h5")
txi = tximport(files, type="kallisto", tx2gene=tx2gene, ignoreAfterBar=TRUE, ignoreTxVersion=TRUE)
write.csv(txi$counts,file="merged_transcripts.csv")

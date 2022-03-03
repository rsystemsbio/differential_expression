#################
# This script takes gene counts from merged Kallisto counts, cleans the noise and low count genes, 
# acquires Differentially expressed genes, and provides a pdfs of PCA visualization
genecount_file = ''

df <- read.csv(genecount_file, header = TRUE, stringsAsFactors = FALSE)
df.new   <- as.data.frame(df)
df.new2  <- df.new
#Optional to remove poorly preforming sample (14)
# df.new2 <- cbind(genes, df.new[-c(14)])
# rownames(df.new2) <- df.new2$genes
rownames(df.new2)  <- df$X
df.new2 = subset(df.new2, select = -c(1) )
colnames(df.new2)  <- paste0("S", 1:29)

## using apply() to iteratively go through each row (as defined by 2nd input; 1=row, 2=col) and sum each column
col.sum  <- apply(df.new2, 2, sum)

df.new3  <- sapply(1:length(col.sum), function(x) df.new2[,x]/col.sum[x]*1e6)
rownames(df.new3)  <- df$X
colnames(df.new3)  <- paste0("S", 1:ncol(df.new3))
head(df.new3)

## using apply() on a per row basis (1 = row in apply()) to determine fraction of elements per row that == zero
zerofrac  <- apply(df.new2, 1, function(x) length(which(x == 0))/ncol(df.new3))

## determine indices of zerofrac (which == row indices) that have a fractional amount of zeros >= 0.95
rows2rm   <- which(zerofrac >= 0.95)

## remove rows that have zero fraction >= 0.95 using the rows2rm (agian, which == indices of those rows that do NOT meet zero threshold)
df.new2  <- df.new2[-rows2rm, ]
df.new3.log  <- log1p(df.new2)
#df.new3.log.zscore1 <- scale(df.new3.log)

## z scores of log1p values using apply function
df.new3.log.zscore  <- apply(df.new3.log, 1, function(x) (x - mean(x))/sd(x))

## rank order genes based on variance
gene.var  <- apply(df.new3.log, 1, var)
var.idx   <- order(gene.var, decreasing = TRUE)
ntop      <- 500

KO.cols  <- c(rep("grey77", 3), rep("green", 3), rep("orange", 3), rep("red", 3), rep("red4", 3), 
              rep("cyan", 3), rep("dodgerblue", 3), rep("blue", 3), rep("magenta", 3), rep("black", 3))

KO.cols.2  <- KO.cols[-14]
meta.data = read.delim(file="", row.names = 1)
meta.data = meta.data[meta.data[-(14), ]]
                      
write.csv(df.new3.log.zscore,file="variance_df_apr9.csv")
cols.to.keep<-which(colnames(df.new3.log.zscore) %in% pc_union)
df.new3.log.zscore <- df.new3.log.zscore[,cols.to.keep]
pdf("CRISPR_KO_Kallisto_TOP_PCs.pdf")
for(i in c(500, 1000, 1500, ncol(df.new3.log.zscore))){
  ntop      <- i
  
  ## PCA minus sample 14
  pca.res     <- prcomp(df.new3.log.zscore[, var.idx[1:ntop]], center = FALSE, scale. = FALSE)
  
  PC_var_exp  <- pca.res$sdev^2/(sum(pca.res$sdev^2))  # this is the eigenvalue of each PC - this indicates how much variation in 
  dim(pca.res$x)
  xlims  <- c(min(pca.res$x[, 1]) -5, max(pca.res$x[, 1]) + 40)
  ylims  <- c(min(pca.res$x[, 2]) -5, max(pca.res$x[, 2]) + 5)
  
  
  plot(pca.res$x[, 1:2], xlim = xlims, ylim = ylims, 
       xlab = paste0("PC1: ", signif(PC_var_exp[1]*100, 3),"% variance"),
       ylab = paste0("PC2: ", signif(PC_var_exp[2]*100, 3),"% variance"),
       # pch = 21, bg = KO.cols, cex = 1.25,
       pch = 21, bg = KO.cols.2, cex = 1.25,
       main = paste0("CRISPR KO (Kallisto-aligned data, DESeq2 norm prior - top Vars ", ntop, " var genes"))
  text(pca.res$x[, 1:2], rownames(df.new3.log.zscore)[-14], cex = 0.5, pos = 4, offset = 0.5)
  #text(pca.res$x[, 1:2], rownames(df.new3.log.zscore)[-14], cex = 0.5, pos = 4, offset = 0.5)
  # text(pca.res$x[, 1:2], sampleTable_v2$sampleName[-14], cex = 0.5, pos = 4, offset = 0.5)
  legend("topright", legend = unique(meta.data$condition), pch = 21, pt.bg = unique(KO.cols), pt.cex = 1.1, cex = 0.5)
}
dev.off()
pca.res     <- prcomp(df.new3.log.zscore[, var.idx[1:40470]], center = FALSE, scale. = FALSE)

pca.res  <- prcomp(t(df.new3.log), center=FALSE, scale.=FALSE)
pc_union = NULL

for(i in 1:16 ) {
  pc_col <- sort(abs(pca.res$rotation[, i]))
  pc_keep <- pc_col[22928:length(pc_col)]
  pc_union <- c(pc_union, names(pc_keep))
}
pc_union <- unique(pc_union)
length(pc_union)

write.csv(pc_union,file="top_100_var_genes.csv")
meta.data = read.delim(file="", row.names = 1)

df.try <- round(df.new3.log)

library(DESeq2)
dds <- DESeqDataSetFromMatrix(countData=round(df.new2), 
                              colData=meta.data, design= ~ condition)

cts <- counts(dds)
geoMeans <- apply(cts, 1, function(row) if (all(row == 0)) 0 else exp(mean(log(row[row != 0]))))
dds <- estimateSizeFactors(dds, geoMeans=geoMeans)
sizeFactors(dds)

colData(dds)$condition<-relevel(dds$condition,ref="")
Uleveldds<-DESeq(dds)
U_B_KO<-results(Uleveldds,contrast=c("condition","",""))
write.csv(U_B_KO,file=".csv")

U_E_KO<-results(Uleveldds,contrast=c("condition","",""))
write.csv(U_E_KO,file=".csv")

U_Ep_KO<-results(Uleveldds,contrast=c("condition","",""))
write.csv(U_Ep_KO,file=".csv")

colData(dds)$condition<-relevel(dds$condition,ref="U_GFP_ffluc")
U_GFP_ffluc<-DESeq(dds)
U_I_KO_GFP_ffluc<-results(U_GFP_ffluc,contrast=c("condition","",""))
write.csv(U_I_KO_GFP_ffluc,file=".csv")

"""
This section removes noisy genes (based on noise threshold) found in either control or treated samples
from the DEG list from DESEq2.
"""
normalized_counts <- counts(dds, normalized=TRUE)

allseqs <- round(df.new2, 0) + 1

#Prepping sample groups
U <- allseqs[,1:3]
U_GFP_ffluc <- allseqs[,4:6]
U_BH73 <- allseqs[,7:14]
U_E <- allseqs[,15:23]
U_I_KO_GFP_ffluc <- allseqs[,24:26]
U_Ep <- allseqs[,27:29]

####################################
#For B
####################################
B_file = '/Users/rachel/U_B_KOf_reverse.csv'
B_deseq2 <- read.csv(B_file, header = TRUE, stringsAsFactors = FALSE)
rownames(B_deseq2) <- B_deseq2$X
B_deseq2 <- B_deseq2[,-1]
B_deseq2 %>% filter(!is.na(padj))
B_deseq2 <- filter(B_deseq2, padj < .05)
B_deseq2_1 <- filter(B_deseq2, log2FoldChange < -1)
B_deseq2__1 <- filter(B_deseq2, log2FoldChange > 1)
B_deseq2 <- rbind(B_deseq2_1, B_deseq2__1)

gene_means     <- rep(NA, nrow(normalized_counts))
gene_sd        <- rep(NA, nrow(normalized_counts))
U_noise     <- rep(NA, nrow(normalized_counts))
names(U_noise) <- rownames(normalized_counts)
U_B_noise     <- rep(NA, nrow(normalized_counts))
names(U_B_noise) <- rownames(normalized_counts)

#BH73 vs U
for(j in 1:nrow(allseqs)){
  gene_means[j]       <- mean(U[j,])
  gene_sd[j]          <- sd(U[j,])
  U_noise[j] <- (gene_sd[j]/sqrt(ncol(U))/(gene_means[j]))
  gene_means[j]       <- mean(U_BH73[j,])
  gene_sd[j]          <- sd(U_BH73[j,])
  U_B_noise[j] <- (gene_sd[j]/sqrt(ncol(U_BH73))/(gene_means[j]))
}
head(U_noise)
head(U)
U_B_comp <- cbind(U_B_noise, U_noise)
U_B_comp_tf <- apply(U_B_comp, 1, function(x) max(x) <= .5)

U_B_idx <- which(U_B_comp_tf == TRUE)
x    <- names(U_B_idx)
#length(which(names(U_B_idx) %in% rownames(B_deseq2)))

#class(names(U_B_idx))
#class()
rownames(B_deseq2)
B_deseq2_clean_names  <- intersect(x, rownames(B_deseq2))
rows.to.keep<-which(rownames(df.new2) %in% pc_union)
B_deseq2_clean <- B_deseq2[B_deseq2_clean_names,]
length(which(B_deseq2_clean$log2FoldChange >1))
write.csv(B_deseq2_clean, "B_deseq2_clean11_6f_reverse.csv")

####################################
#For E
####################################
E_file = '/Users/rachel/Downloads/KNORM_FILT_9_22/U_E_KOf_reverse.csv'
E_deseq2 <- read.csv(E_file, header = TRUE, stringsAsFactors = FALSE)
E_deseq2
rownames(E_deseq2) <- E_deseq2$X
E_deseq2 <- E_deseq2[,-1]
E_deseq2 %>% filter(!is.na(padj))
E_deseq2 <- filter(E_deseq2, padj < .05)
E_deseq2_1 <- filter(E_deseq2, log2FoldChange < -1)
E_deseq2__1 <- filter(E_deseq2, log2FoldChange > 1)
E_deseq2 <- rbind(E_deseq2_1, E_deseq2__1)

gene_means     <- rep(NA, nrow(normalized_counts))
gene_sd        <- rep(NA, nrow(normalized_counts))
U_noise     <- rep(NA, nrow(normalized_counts))
names(U_noise) <- rownames(normalized_counts)
U_E_noise     <- rep(NA, nrow(normalized_counts))
names(U_E_noise) <- rownames(normalized_counts)

#E vs U
for(j in 1:nrow(normalized_counts)){
  gene_means[j]       <- mean(U[j,])
  gene_sd[j]          <- sd(U[j,])
  U_noise[j] <- (gene_sd[j]/sqrt(ncol(U))/(gene_means[j]))
  gene_means[j]       <- mean(U_E[j,])
  gene_sd[j]          <- sd(U_E[j,])
  U_E_noise[j] <- (gene_sd[j]/sqrt(ncol(U_E))/(gene_means[j]))
}
head(U_noise)
head(U)
U_E_comp <- cbind(U_E_noise, U_noise)
U_E_comp_tf <- apply(U_E_comp, 1, function(x) max(x) <= .5)

U_E_idx <- which(U_E_comp_tf == TRUE)
x    <- names(U_E_idx)

E_deseq2_filtered  <- intersect(x, rownames(E_deseq2))
E_deseq2_clean <- E_deseq2[E_deseq2_filtered,]
write.csv(E_deseq2_clean, "E_deseq2_cleanf_reverse.csv")
####################################
#For U_I_KO_GFP_ffluc
####################################
I_file = '/Users/rachel/Downloads/KNORM_FILT_9_22/U_I_KO_GFP_fflucf_reverse.csv'
I_deseq2 <- read.csv(I_file, header = TRUE, stringsAsFactors = FALSE)
rownames(I_deseq2) <- I_deseq2$X
I_deseq2 <- I_deseq2[,-1]
I_deseq2 %>% filter(!is.na(padj))
I_deseq2 <- filter(I_deseq2, padj < .05)
I_deseq2_1 <- filter(I_deseq2, log2FoldChange < -1)
I_deseq2__1 <- filter(I_deseq2, log2FoldChange > 1)
I_deseq2 <- rbind(I_deseq2_1, I_deseq2__1)

gene_means     <- rep(NA, nrow(normalized_counts))
gene_sd        <- rep(NA, nrow(normalized_counts))
U_GFP_ffluc_noise    <- rep(NA, nrow(normalized_counts))
names(U_GFP_ffluc_noise) <- rownames(normalized_counts)
U_I_KO_GFP_ffluc_noise     <- rep(NA, nrow(normalized_counts))
names(U_I_KO_GFP_ffluc_noise) <- rownames(normalized_counts)

#I vs U_GFP_Flucc
for(j in 1:nrow(normalized_counts)){
  gene_means[j]       <- mean(U_GFP_ffluc[j,])
  gene_sd[j]          <- sd(U_GFP_ffluc[j,])
  U_GFP_ffluc_noise[j] <- (gene_sd[j]/sqrt(ncol(U_GFP_ffluc))/(gene_means[j]))
  
}
dim(U_I_KO_GFP_ffluc)
U_I_KO_GFP_ffluc
dim(U_I_KO_GFP_ffluc)
for(j in 1:nrow(normalized_counts)){
  gene_means[j]       <- mean(U_I_KO_GFP_ffluc[j,])
  gene_sd[j]          <- sd(U_I_KO_GFP_ffluc[j,])
  U_I_KO_GFP_ffluc_noise[j] <- (gene_sd[j]/sqrt(ncol(U_I_KO_GFP_ffluc))/(gene_means[j]))
}
head(U_noise)
head(U)
U_I_comp <- cbind(U_I_KO_GFP_ffluc_noise, U_GFP_ffluc_noise)
U_I_comp_tf <- apply(U_I_comp, 1, function(x) max(x) <= .5)

U_I_idx <- which(U_I_comp_tf == TRUE)
x    <- names(U_I_idx)

I_deseq2_filtered  <- intersect(x, rownames(I_deseq2))
I_deseq2_clean <- I_deseq2[I_deseq2_filtered,]
write.csv(I_deseq2_clean, "I_deseq2_clean_reverse.csv")

####################################
#For Ep
####################################
Ep_file = '/Users/rachel/Downloads/KNORM_FILT_9_22/U_E_KOf_reverse.csv'
Ep_deseq2 <- read.csv(Ep_file, header = TRUE, stringsAsFactors = FALSE)
rownames(Ep_deseq2) <- Ep_deseq2$X
Ep_deseq2 <- Ep_deseq2[,-1]
Ep_deseq2 %>% filter(!is.na(padj))
Ep_deseq2 <- filter(Ep_deseq2, padj < .05)
Ep_deseq2_1 <- filter(Ep_deseq2, log2FoldChange < -1)
Ep_deseq2__1 <- filter(Ep_deseq2, log2FoldChange > 1)
Ep_deseq2 <- rbind(Ep_deseq2_1, Ep_deseq2__1)

gene_means     <- rep(NA, nrow(normalized_counts))
gene_sd        <- rep(NA, nrow(normalized_counts))
U_noise     <- rep(NA, nrow(normalized_counts))
names(U_noise) <- rownames(normalized_counts)
U_Ep_noise     <- rep(NA, nrow(normalized_counts))
names(U_Ep_noise) <- rownames(normalized_counts)

#Ep vs U
for(j in 1:nrow(normalized_counts)){
  gene_means[j]       <- mean(U[j,])
  gene_sd[j]          <- sd(U[j,])
  U_noise[j] <- (gene_sd[j]/sqrt(ncol(U))/(gene_means[j]))
  gene_means[j]       <- mean(U_Ep[j,])
  gene_sd[j]          <- sd(U_Ep[j,])
  U_Ep_noise[j] <- (gene_sd[j]/sqrt(ncol(U_Ep))/(gene_means[j]))
}
head(U_noise)
head(U)
U_Ep_comp <- cbind(U_Ep_noise, U_noise)
U_Ep_comp_tf <- apply(U_Ep_comp, 1, function(x) max(x) <= .5)

U_Ep_idx <- which(U_Ep_comp_tf == TRUE)
x    <- names(U_Ep_idx)

Ep_deseq2_filtered  <- intersect(x, rownames(Ep_deseq2))
Ep_deseq2_clean <- Ep_deseq2[Ep_deseq2_filtered,]
write.csv(Ep_deseq2_clean, ".csv")


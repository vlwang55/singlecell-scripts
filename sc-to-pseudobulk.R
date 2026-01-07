# This method uses DeSEQ2 as it is more accurate than using the Seurat FindMarkers
# It requires at least 2 replicates in at least one identity you are comparing (i.e. genotype needs 2 WTs and at least 1 NULL)
# Ideally you want at least 2 replicates in all of the identities you are comparing, but it's single cell data so that doesn't always happpen

##### MAKE PSEUDOBULK COUNT MATRIX AND METADATA TABLE #####
# Load libraries
library(Seurat)
library(Signac)
library(dplyr)
library(ggplot2)
library(tidyverse)
library(Matrix)

DefaultAssay(int) <- 'RNA'

# OPTIONAL For integrated samples
int[["RNA"]] <- JoinLayers(int[["RNA"]])
Idents(int) <- int$orig.ident  # or whatever identity you're interested in
int_list <- SplitObject(int, split.by = "orig.ident")

# Make pseudobulk counts matrix for each sample (can replace the list with just the Seurat object of one sample if only one)
pseudobulk_counts <- lapply(int_list, function(x) {
  counts <- GetAssayData(x, layer = "counts")
  cell <- Idents(x)
  rowsum(as.matrix(t(counts)), group = cell)
})

# Create a combined matrix of all samples together
combined_mat <- do.call(cbind, lapply(seq_along(pseudobulk_counts), function(i) {
  mat <- t(pseudobulk_counts[[i]])  # now genes are rows
  sample_id <- names(paste(merge_list))[i]  # get sample name # rename cols
  return(mat)
}))

# If the above doesn't work, then use the code below although it's less neat
combined_mat <- cbind(t(pseudobulk_counts[[1]]),t(pseudobulk_counts[[2]]),t(pseudobulk_counts[[3]]),t(pseudobulk_counts[[4]]),t(pseudobulk_counts[[5]]),t(pseudobulk_counts[[6]]))

# Check result
dim(combined_mat)  # genes x (cluster-sample combinations)

# Make metadata table with the categories you're interested in comparing
cluster_sample <- colnames(combined_mat)
metadata_table <- data.frame(
 genotype = c("WT", "CORR", "NULL", "NULL"),
  batch = c('2', '2', '2', '1')
  )
rownames(metadata_table) <- cluster_sample

# Order your identity levels based on what you want the DESeq to compare - the first one listed will always be compared to the others
metadata_table$genotype <- factor(metadata_table$funcgeno, levels = c("NULL", "CORR", "WT"))


##### RUN DESEQ2 #####
# Load library
library(DESeq2)

dds <- DESeqDataSetFromMatrix(
  countData = combined_mat,
  colData = metadata_table,
  design = ~genotype) # the metadata category you want to compare

# OPTIONAL if you just want the normalized matrix of counts for heatmaps then use this 
rld <- rlog(dds, blind=TRUE)
normalized_matrix <- assay(rld)

# Again this will only work if you have at least 2 replicates for one of the identities in genotype
dds <- DESeq(dds)
vsd <- vst(dds, blind = TRUE)

# OPTIONAL For batch correction if you have multiple batches
batch_corrected <- removeBatchEffect(assay(vsd), batch = vsd$batch)
assay(vsd) <- batch_corrected

# Make PCA plot
DESeq2::plotPCA(vsd, intgroup = "geno")


##### DEG #####
# lists the coefficients
resultsNames(vsd) 
res <- results(dds, name="genotype_NULL_vs_WT")
write.table(res, file = "PseudobulkBatchS4allParseallGenes_NullNotchinhvsNull.csv", sep=",", quote =F, col.names=NA)


##### HEATMAPS #####
library(pheatmap)
library(RColorBrewer)

matrix <- assay(vsd)

# Ensure genes are in the matrix or else it will error 
genelist %in% rownames(matrix)

# Heat map of gene expressions with row clustering
pheatmap(matrix[genelist,],
         cluster_rows=R,
         show_rownames=TRUE,
         cluster_cols=F,
         scale = 'none',
         color=colorRampPalette(c("blue ","white","red"))(100))


##### TO RENORMALIZE SAMPLES VIA TPM (useful if comparing to other BulkRNAseq data) #####
# Make a tpm function
tpm <- function(counts, gene_length) {
  # Normalize counts by gene length (per kilobase)
  rate <- counts / (gene_length / 1000)
  
  # Calculate the "per million" scaling factor
  denom <- sum(rate) / 1e6
  
  # Calculate TPM
  xtpm <- rate / denom
  
  return(xtpm)
}

# Assuming 'counts_df' is a matrix/data frame of counts (genes as rows, samples as columns)
# and 'gene_lengths' is a matching vector/data frame of gene lengths.

library(dplyr)

TPM_matrix <- apply(X = combined_mat, MARGIN = 2, FUN = function(x) {
tpm(x, gene_lengths$length)
}) %>% as.data.frame()

# If you want to combine with a different bulk RNA dataset, load your TPM values for the other dataset(s)
abundance2 <- read.table("H1S0_TPM.csv", sep = ",")
rownames(abundance2) <- abundance2$V1
dim(abundance2)
dim(tpm.mat)
abundance3 <- abundance2[rownames(tpm.mat),]
dim(abundance3)
head(tpm.mat)

tpmfinal <- cbind(tpm.mat, 'obj' = as.numeric(abundance3$V2))
head(tpmfinal)

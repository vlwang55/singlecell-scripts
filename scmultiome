##### This script is a modified version of https://github.com/quadbio/scMultiome_analysis_vignette/blob/main/Tutorial.md
##### Also contains cell cycle scoring and motif analysis from https://satijalab.org/seurat/archive/v4.3/weighted_nearest_neighbor_analysis

# Installation
BiocManager::install("TFBSTools", type = "source", force = TRUE)
install.packages("TFBSTools")

# Load libraries
library(Seurat)
library(EnsDb.Hsapiens.v86)
library(Signac)
library(dplyr)
library(ggplot2)
library(GenomicRanges)
library(tidyverse)
library(GenomeInfoDb)
library(BSgenome.Hsapiens.UCSC.hg38)
library(remotes)
library(TFBSTools)


##### READING IN YOUR FILES #####
# 10x hdf5 file 
inputdata.10x <- Read10X_h5("filtered_feature_bc_matrix.h5")
fragpath <- "filelocation/atac_fragments.tsv.gz"

# Extract RNA and ATAC data
rna_counts <- inputdata.10x$`Gene Expression`
atac_counts <- inputdata.10x$Peaks

# Create Seurat object and calculate mitochondrial gene percentage for each cell - for mouse it is ^mt-
# DefaultAssay is how you set what assay you're working in. You can check your assay just by doing DefaultAssay(obj)
obj <- CreateSeuratObject(counts = rna_counts, project = "ProjName")
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
DefaultAssay(obj) <- "RNA"

# Create ATAC assay and add it to the object - grange associates the chromosome locations with gene names
grange.counts <- StringToGRanges(rownames(atac_counts), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac_counts <- atac_counts[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"

obj[["ATAC"]] <- CreateChromatinAssay(
  counts = atac_counts,
  sep = c(":", "-"),
  genome = 'hg38',
  fragments = fragpath,
  annotation = annotations
)

# Clean up your workspace - remove gets rid of objects and gc() refreshes your working data
remove(inputdata.10x, rna_counts, atac_counts, grange.counts, grange.use)
gc()

##### PEAK CALLING USING MACS3 ##### 
# For this section you will need to download MACS3 and a genomic blacklist region list ahead of time
# MACS3 identifies peaks - regions with significant counts
DefaultAssay(obj) <- "ATAC"
peaks <- CallPeaks(obj,
                   macs2.path = "macs3filelocation/macs3",
                   effective.genome.size= 2.7e9)

# Remove peaks on nonstandard chromosomes and in genomic blacklist regions
setwd("genomicblacklistfilelocation/")
blacklist_hg38 <- read.table("blacklist.v2.bed", header = FALSE, sep = "\t")
colnames(blacklist_hg38) <- c("Chr","Start", "End", "Quality?")
blacklist_hg38 <- makeGRangesFromDataFrame(blacklist_hg38)
peaks <- keepStandardChromosomes(peaks, pruning.mode = "coarse")
peaks <- subsetByOverlaps(x = peaks, ranges = blacklist_hg38, invert = TRUE)
remove(blacklist_hg38)

# OPTIONAL  ### How to rewrite fragment file locations if you switch servers ###
# Your fragment file location or fragpath needs to be updated if you move between servers
newfragobj <- CreateFragmentObject(path = "newfilelocation/atac_fragments.tsv.gz", cells = colnames(obj), validate.fragments = TRUE)
Fragments(obj@assays$ATAC) <- newfragobj

##### QUANTIFY COUNTS PER PEAK #####
macs3_counts <- FeatureMatrix(
  fragments = Fragments(obj),
  features = peaks,
  cells = colnames(obj))

# Create a new assay using the MACS2 peak set and add it to the Seurat object
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
obj[["peaks"]] <- CreateChromatinAssay(
  counts = macs3_counts,
  fragments = Fragments(obj),
  annotations = annotations)

##### QUALITY CONTROL FOR READS #####
# Get rid of low quality cells
# These numbers are guidelines for our iPSC diff pancreatic cells - other cell types will vary
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA","nCount_peaks","percent.mt"),
        log = TRUE, pt.size = 0, group.by= "orig.ident") + NoLegend()

obj <- subset(
  x = obj,
  subset = nFeature_RNA > 500 & 
    nFeature_RNA < 7000 & 
    nCount_RNA > 1500 & 
    nCount_RNA < 50000 & 
    percent.mt < 15 &
    nCount_peaks > 5000 &
    nCount_Peaks < 100000
)

# Double check that your cells have been filtered
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA","nCount_peaks","percent.mt"),
        log = TRUE, pt.size = 0, group.by= "orig.ident") + NoLegend()

# Save your file!! Do this often
saveRDS(obj, file ="MyObj.rds")


##### RNA ANALYSIS #####
# This method uses SCTransform which is helpful because it preserves the original RNA counts in the RNA assay
# This is very useful if you want to later subset a portion of the object and renormalize it
DefaultAssay(obj) <- "RNA"
obj <- SCTransform(obj, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap', reduction.key = 'rnaUMAP_')

##### ATAC ANALYSIS #####
# We exclude the first dimension as this is typically correlated with sequencing depth
DefaultAssay(obj) <- "peaks"
obj <- RunTFIDF(obj)
obj <- FindTopFeatures(obj, min.cutoff = 'q0')
obj <- RunSVD(obj)
obj <- RunUMAP(obj, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

##### WNN ANALYSIS AND UMAP GENERATION #####
# Weighted nearest neighbor analysis integrates both ATAC and RNA data and weights it for each cell for more accurate clustering
obj <- FindMultiModalNeighbors(obj, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
obj <- RunUMAP(obj, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
obj <- FindClusters(obj, graph.name = "wsnn", algorithm = 3, verbose = FALSE)

##### VISUALIZE UMAPS #####
p1 <- DimPlot(obj, reduction = "wnn.umap", label = TRUE, repel = TRUE) + ggtitle("WNN TITLE")
p2 <- DimPlot(obj, reduction = "umap.atac", label = TRUE, repel = TRUE) + ggtitle("ATAC")
p3 <- DimPlot(obj, reduction = "umap.rna", label = TRUE, repel = TRUE) + ggtitle("RNA")
p1 + p2 + p3 & NoLegend() & theme(plot.title = element_text(hjust = 0.5))

##############################

##### ASSIGN CLUSTER IDENTITIES #####
# You don't have to use all of these methods - mix and match depending on your celltypes

### Show RNA expression of gene markers to help assign cluster identities
# In my experience this is more important than looking at the top DEGs in each cluster
DefaultAssay(obj) <- "RNA"
FeaturePlot(obj, features = c('PDX1', 'NKX6-1', 'NEUROD1', 'COL3A1'), label = TRUE, pt.size=0.5)

### Calculate top DEGs for each cluster
ClusterBiomarkers <- data.frame(matrix(nrow=500))
obj <- PrepSCTFindMarkers(obj)
# Change based on the number of clusters you have
for (x in 0:6) {
  cluster.markers <- FindMarkers(obj, ident.1 = (x), min.pct = 0.25)
  cluster.markers_p05 <- subset(cluster.markers, p_val < 0.05 & avg_log2FC > 0)
  cluster.markers_p05_desc <- arrange(cluster.markers_p05, desc(avg_log2FC))
  ClusterBiomarkers<-cbind(ClusterBiomarkers, rownames(cluster.markers_p05_desc[1:500,]))
}

ClusterBiomarkers <- ClusterBiomarkers[,-1]
colnames(ClusterBiomarkers) <- 0:6
write.table(ClusterBiomarkers, file = "ClusterDEGs.csv", sep=",", quote =T, row.names=F)
# Now can input it into http://bio-bigdata.hrbmu.edu.cn/CellMarker/CellMarker_annotation.jsp

### Find all markers distinguishing cluster 5 from clusters 0 and 3
cluster5.markers <- FindMarkers(obj, ident.1 = '5', ident.2 = c('0', '3'), min.pct = 0.25)
cluster5.markers_desc <- arrange(cluster5.markers, desc(avg_log2FC))
head(cluster5.markers_desc, n = 50)

### Find markers for every cluster compared to all remaining cells, report only the positive ones
cluster2.markers <- FindMarkers(obj, ident.1 = 2, logfc.threshold = 0.25, test.use = "roc", only.pos = TRUE)
cluster2.markers_desc <- arrange(cluster2.markers, desc(avg_log2FC))
head(cluster2.markers_desc, n = 50)

### Perform sub-clustering on cluster 6 to find additional structure
obj <- FindSubCluster(obj, cluster = "6", graph.name = 'snn', resolution = 0.1)
Idents(obj) <- "sub.cluster"
DimPlot(obj, label=TRUE)

# Add cluster annotations and save under meta.data as celltypes
obj <- RenameIdents(obj, '0' = 'PancProg', '1' = 'PreEndocrine', '2' = 'Endothelial', '3' = 'Acinar')
obj@meta.data$celltypes <- Idents(obj)

##### ATAC DATA VISUALIZATION ##### 
# To make the visualization easier, can subset cell clusters
celltype.names <- levels(seu)
maincluster.names <- grep("GATA6|proliferative|Neuron|peri", celltype.names,value = TRUE)
maincluster <- subset(seu_rename, idents = maincluster.names)
CoveragePlot(CCPP, region = 'chr8-9903142-9903717', assay = 'ATAC',expression.assay = 'integrated', peaks = TRUE)

##############################

##### OTHER ######
# Number of cells per cluster
table(obj@active.ident, obj@meta.data$orig.ident)

# Heatmap and violin plot for RNA expression
VlnPlot(obj, features = c('PDX1'), group.by = 'orig.ident')
VlnPlot(obj, group.by = "orig.ident", features = c('PDX1'), stack = TRUE)
DoHeatmap(obj, features = top50DEG, group.by= "orig.ident")

# Cell cycle scoring
library(Seurat)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
obj <- CellCycleScoring(obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

RidgePlot(obj, features = c("PCNA", "TOP2A", "MCM6", "MKI67"), ncol = 2)

obj <- RunPCA(obj, features = c(s.genes, g2m.genes))
DimPlot(obj, label= TRUE, reduction = "pca")

# Cell cycle scoring for mice
library(dplyr)
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

homologs = read.csv("http://www.informatics.jax.org/downloads/reports/HOM_MouseHumanSequence.rpt",sep="\t")
ms2hs <- function(gene_list){
  output = c()
  for(gene in gene_list){
    class_key = (homologs %>% filter(Symbol == gene & Common.Organism.Name=="human"))[['DB.Class.Key']]
    if(!identical(class_key, integer(0)) ){
      human_genes = (homologs %>% filter(DB.Class.Key == class_key & Common.Organism.Name=="mouse, laboratory"))[,"Symbol"]
      for(human_gene in human_genes){
        output = append(output,human_gene)
      }
    }
  }
  return (output)
}

s.genes <- ms2hs(s.genes)
g2m.genes <- ms2hs(g2m.genes)

mouseobj <- CellCycleScoring(mouseobj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

##### FINDING MOTIFS #####
library(chromVAR)
library(JASPAR2020)
library(TFBSTools)
library(motifmatchr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Scan the DNA sequence of each peak for the presence of each motif, and create a Motif object
DefaultAssay(obj) <- "peaks"
pwm_set <- getMatrixSet(x = JASPAR2020, opts = list(species = 9606, all_versions = FALSE))
motif.matrix <- CreateMotifMatrix(features = granges(obj), pwm = pwm_set, genome = "hg38", use.counts = FALSE)
motif.object <- CreateMotifObject(data = motif.matrix, pwm = pwm_set)
obj <- SetAssayData(obj, assay = 'peaks', slot = 'motifs', new.data = motif.object)

# Note that this step can take 30-60 minutes 
obj <- RunChromVAR(
  object = obj,
  genome = BSgenome.Hsapiens.UCSC.hg38
)
saveRDS(obj,"objwithmotif.rds")

# Finds motif ID for gene and plots expression of motif alongside RNA expression of gene
motif.name <- ConvertMotifID(obj, name = 'THRB')
gene_plot <- FeaturePlot(obj, features = "sct_THRB", reduction = 'umap')
motif_plot <- FeaturePlot(obj, features = motif.name, min.cutoff = 0, cols = c("lightgrey", "darkred"), reduction = 'umap')
gene_plot | motif_plot
gene_plot

# Differential expression for the motif across all cell types
remotes::install_github("immunogenomics/presto")
markers_rna <- presto:::wilcoxauc.Seurat(X = obj, group_by = 'celltype', assay = 'data', seurat_assay = 'RNA')
markers_motifs <- presto:::wilcoxauc.Seurat(X = obj,, group_by = 'celltype', assay = 'data', seurat_assay = 'chromvar')
motif.names <- markers_motifs$feature
colnames(markers_rna) <- paste0("RNA.", colnames(markers_rna))
colnames(markers_motifs) <- paste0("motif.", colnames(markers_motifs))
markers_rna$gene <- markers_rna$RNA.feature
markers_motifs$gene <- ConvertMotifID(obj, id = motif.names)

# a simple function to implement the procedure above
topTFs <- function(celltype, padj.cutoff = 1e-2) {
  ctmarkers_rna <- dplyr::filter(
    markers_rna, RNA.group == celltype, RNA.padj < padj.cutoff, RNA.logFC > 0) %>% 
    arrange(-RNA.auc)
  ctmarkers_motif <- dplyr::filter(
    markers_motifs, motif.group == celltype, motif.padj < padj.cutoff, motif.logFC > 0) %>% 
    arrange(-motif.auc)
  top_tfs <- inner_join(
    x = ctmarkers_rna[, c(2, 11, 6, 7)], 
    y = ctmarkers_motif[, c(2, 1, 11, 6, 7)], by = "gene"
  )
  top_tfs$avg_auc <- (top_tfs$RNA.auc + top_tfs$motif.auc) / 2
  top_tfs <- arrange(top_tfs, -avg_auc)
  return(top_tfs)
}

# Identify top markers in any of the cell types and visualize
top3TF <- data.frame(matrix(ncol = 9, nrow = 0))  
newidents <- c('PancProg','ProlifGATA6','GATA6','Mesenchyme', 'PreEndo', 'EMTlike','Stellate')

for (x in newidents) {
  topTFtemp <- topTFs(x)
  top3TF <- rbind(top3TF, topTFtemp[1:3, ])  # Append rows to top3TF
}

write.table(top3TF, file = "top3TF.csv", sep=",", quote =T, row.names=F)

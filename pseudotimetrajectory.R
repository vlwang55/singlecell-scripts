# This script is made up of modified versions of the Monocle3 and Slingshot vignettes

##### MONOCLE #####
# Load libraries
library(Signac)
library(monocle3)
library(SeuratWrappers)

# Convert Seurat object to cds and process
DefaultAssay(obj) <- "SCT"
obj.cds <- as.cell_data_set(obj)
obj.cds <- preprocess_cds(obj.cds, method = 'PCA')
obj.cds <- reduce_dimension(obj.cds, reduction_method = 'UMAP')
obj.cds <- cluster_cells(cds = obj.cds, reduction = "UMAP")
obj.cds <- learn_graph(obj.cds, use_partition = TRUE)

# Order cells
obj.cds <- order_cells(NN3.cds, reduction_method = "UMAP")

obj.cds<- choose_graph_segments(NN3.cds, reduction_method = "UMAP")


# Plot trajectories colored by pseudotime
plot_cells(
  cds = obj.cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE,
)

# Note I have not tested this section thoroughly
n # a helper function to identify the root principal points:
get_earliest_principal_node <- function(obj.cds, time_bin="130-170"){
  cell_ids <- which(colData(NN3.cds)[, "embryo.time.bin"] == time_bin)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

obj.cds <- order_cells(obj.cds, root_pr_nodes=get_earliest_principal_node(obj.cds))


##### Slingshot #####
# Load libraries
library(slingshot)
library(grDevices)
library(RColorBrewer)

# First convert Seurat object into a Single Cell Experiment
sce <- as.SingleCellExperiment(merge, assay = "RNA")

# Run Slingshot - can do on PCA instead of UMAP as well
sce <- slingshot(sce, clusterLabels = 'celltypes', reducedDim = 'UMAP')

# Set color palette
colors <- colorRampPalette(brewer.pal(11,'Spectral')[-6])(100)
plotcol <- colors[cut(sce$slingPseudotime_1, breaks=100)]

# Visualize
plot(reducedDims(sce)$UMAP, col = plotcol, pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, col='black')
plot(reducedDims(sce)$UMAP, col = brewer.pal(9,'Set1')[sce$celltypes], pch=16, asp = 1)
lines(SlingshotDataSet(sce), lwd=2, type = 'lineages', col = 'black')

# This script is a modified version of https://satijalab.org/seurat/articles/integration_introduction
# Can be used both for RNA/ ATAC or just RNA

##### SECTION 1: RNA ONLY #####
library(glmGamPoi)
library(Seurat)
library(SeuratData)

# OPTIONAL might need to increase working data
options(future.globals.maxSize= 891289600)

# Start with your Seurat objects that have already been filtered for quality
mylist <- c(obj1, obj2)

# It is very important that we use SCTransform for normalization 
# SCTransform preserves the initial RNA raw counts in the RNA assay
mylist <- lapply(X = mylist, FUN = function(x) {
  x <- SCTransform(x)})

# Select features that are repeatedly variable across datasets for integration and then integrate
features <- SelectIntegrationFeatures(object.list = mylist)
anchors <- FindIntegrationAnchors(object.list = mylist, anchor.features = features)
int <- IntegrateData(anchorset = anchors)

# You should now have an 'integrated' assay and can run downstream processing as normal
DefaultAssay(int) <- "integrated"
int <- ScaleData(int, verbose = F)
int <- RunPCA(int, npcs = 50, verbose = F)
int <- RunUMAP(int, reduction = "pca", dims = 1:20, verbose = F)
DimPlot(int)
saveRDS(int, "RNA_int.rds")

##### SECTION 2: RNA AND ATAC #####
# Again, make sure you finished peak calling and filtering your objects
# Do as above but now also do the following:

DefaultAssay(int) <- "peaks"
newlist <- SplitObject(int, split.by = "orig.ident")

# You will have to run normalization and variable feature finding on both your split list and your merged object
int <- FindTopFeatures(int, min.cutoff = 10)
int <- RunTFIDF(int)
int <- RunSVD(int)

newlist <- lapply(X = newlist, FUN = function(x) {
  x <- FindTopFeatures(x, min.cutoff = 10)
  x <- RunTFIDF(x)
  x <- RunSVD(x)
  })

# find integration anchors
integration.anchors <- FindIntegrationAnchors(
  object.list = newlist,
  anchor.features = rownames(int),
  reduction = "rlsi",
  dims = 2:30
)

# integrate LSI embeddings
integrate2 <- IntegrateEmbeddings(
  anchorset = integration.anchors,
  reductions = int[["lsi"]],
  new.reduction.name = "integrated_lsi",
  dims.to.integrate = 1:30
)

int[['integrated_lsi_peaks']] <- CreateDimReducObject(
  Embeddings(integrate2, "integrated_lsi")[colnames(int),], key="INTEGRATEDLSIPEAKS_", assay="peaks"
)

int <- RunUMAP(int,
                 reduction = "integrated_lsi_peaks",
                 dims = 2:30,
                 reduction.name = "umap_seurat_peaks",
                 reduction.key = "UMAPSEURATpeaks_")
# Test it by looking at UMAP, then move on to WNN

# Now can do WNN
DefaultAssay(int) <- "integrated"
int <- FindMultiModalNeighbors(int, reduction.list = list("pca", "integrated_lsi_peaks"), dims.list = list(1:30, 2:50), modality.weight.name = c("RNA.weight","ATAC.weight"),
                                 verbose = TRUE)
int <- RunUMAP(int, nn.name = "weighted.nn", assay = "RNA")
int <- FindClusters(int, graph.name = "umap", resolution = 0.1)
UMAPPlot(merge, label = TRUE) & NoAxes()

saveRDS(int, "RNA_ATAC_int.rds")


##### SECTION 3: RNA WHERE YOU HAVE SUBSETTED CLUSTERS #####
# Now let's say you've generated the RNA integrated object from Section 1 above (or just an RNA integrated object)
# but you want to subset out a cluster and renormalize and reintegrate
# This is tricky because you need to get rid of all your assays except for RNA which has the raw counts thanks to SCTransform

# First subset your identity(s) of interest
sub <- subset(int, idents = '1')

DefaultAssay(sub) <- 'RNA'
sub[['SCT']] <- NULL
sub[['integrated']] <- NULL

# You will be left with raw counts in the RNA assay split by their original samples- you need to combine these into one matrix
# before you can split the matrix by sample- otherwise you will have different numbers of genes
sub[["RNA"]] <- JoinLayers(sub[["RNA"]])
newnewlist <- SplitObject(sub, split.by = "orig.ident")

# Now we can do NormalizeData, although you can still do SCTransform if you'd like
# Mostly I do NormalizeData because it is easier to control nfeatures for variable features
newnewlist <- lapply(X = newnewlist, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000)
})

# Redo integration based on variable feature anchors and downstream processing
features <- SelectIntegrationFeatures(object.list = newnewlist)
anchors <- FindIntegrationAnchors(object.list = newnewlist, anchor.features = features)
merge2 <- IntegrateData(anchorset = anchors)
merge2 <- ScaleData(merge2)
merge2 <- RunPCA(merge2)
merge2 <- RunUMAP(merge2, reduction = "pca", dims = 1:20, verbose = F)
merge2 <- FindNeighbors(merge2, dims = 1:20)
merge2 <- FindClusters(merge2, resolution = 0.1)
DimPlot(merge2)

saveRDS(merge2, 'submerge.rds')


# If you need to subset and renormalize with RNA and ATAC objects, then I think it should work the same but just
# get rid of the peaks assay as well and check to see if the ATAC assay needs to have its layers combined like
# RNA did. Then it is a repeat of Section 2. Good luck!

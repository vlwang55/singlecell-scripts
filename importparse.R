## If you haven't I recommend working through the scMultiome script first

# Load libraries
library(Seurat)
library(dplyr)
library(Matrix)
library(ggplot2)

# First, download your Parse files - if they saved as gz, open the csv files (all_genes.csv and cell_metadata.csv) in Notepad and save them as Normal Text File with a .csv
# For the mtx file, open in Notebook and save as "All types" with .mtx

##### READING IN YOUR FILES #####
mat_path <- "/path/to/file"
mat <- ReadParseBio(mat_path)

# Check to see if empty gene names are present, add name if so
table(rownames(mat) == "")
rownames(mat)[rownames(mat) == ""] <- "unknown"

# Read in cell meta data
cell_meta <- read.csv(paste0(mat_path, "/cell_metadata.csv"), row.names = 1)

# Create object, note that the Parse consultant recommended min.cells and min.features of 0, and the Parse website recommends 100. This is what I found worked.
obj <- CreateSeuratObject(mat, project = "ProjName",min.cells = 1, min.features = 250, names.field = 0, meta.data = cell_meta)

##### QUALITY CONTROL #####
obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3,
        log = TRUE, pt.size = 0, group.by= "orig.ident") + NoLegend()

# Filter cells - exact ranges will vary based on cell type
obj <- subset(
  x = obj,
  subset = nFeature_RNA < 5500 & 
    nFeature_RNA > 500 & 
    nCount_RNA < 50000 & 
    nCount_RNA > 1000 & 
    percent.mt < 10)

saveRDS(obj, 'obj_filtered.rds')

# Continue with integration or normalization as normal

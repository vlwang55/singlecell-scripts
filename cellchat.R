### Modified version of https://github.com/sqjin/CellChat/blob/master/tutorial/CellChat-vignette.Rmd

# Load libraries
library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)

obj <- readRDS("obj.rds")

##### CREATE CELLCHAT OBJECT FROM SEURAT OBJECT #####
DefaultAssay(obj) <- "SCT"
cellchat <- createCellChat(object = obj, group.by = "orig.ident")

### OPTIONAL If cell meta information is not added when creating CellChat object, USERS can also add it later using addMeta, and set the default cell identities using setIdent.
cellchat <- addMeta(cellchat, meta = meta)
cellchat <- setIdent(cellchat, ident.use = "labels") # set "labels" as default cell identity
levels(cellchat@idents) # show factor levels of the cell labels
groupSize <- as.numeric(table(cellchat@idents)) # number of cells in each cell group

# Set database
CellChatDB <- CellChatDB.human # use CellChatDB.mouse if running on mouse data
showDatabaseCategory(CellChatDB)

# Show the structure of the database
dplyr::glimpse(CellChatDB$interaction)

# Use a subset of CellChatDB for cell-cell communication analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
# OR use all CellChatDB for cell-cell communication analysis
CellChatDB.use <- CellChatDB # simply use the default CellChatDB

# set the used database in the object
cellchat@DB <- CellChatDB.use

##### RUNNING CELLCHAT ANALYSIS #####
# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) # do parallel

cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# project gene expression data onto PPI (Optional: when running it, USER should set `raw.use = FALSE` in the function `computeCommunProb()` in order to use the projected data)
# cellchat <- projectData(cellchat, PPI.human)

# Note this step takes a while (~30 min)
cellchat <- computeCommunProb(cellchat)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat@netP$pathways

# OPTIONAL We provide a function subsetCommunication to easily access the inferred cell-cell communications of interest. For example,
#df.net <- subsetCommunication(cellchat) returns a data frame consisting of all the inferred cell-cell communications at the level of ligands/receptors. Set slot.name = "netP" to access the the inferred communications at the level of signaling pathways
#df.net <- subsetCommunication(cellchat, sources.use = c(1,2), targets.use = c(4,5)) gives the inferred cell-cell communications sending from cell groups 1 and 2 to cell groups 4 and 5.
#df.net <- subsetCommunication(cellchat, signaling = c("WNT", "TGFb")) gives the inferred cell-cell communications mediated by signaling WNT and TGFb.

df.net <- subsetCommunication(cellchat, signaling = c("NOTCH"))
df.net

##### VISUALIZATION #####
# We can calculate the aggregated cell-cell communication network by counting the number of links or summarizing the communication probability. USER can also calculate the aggregated network among a subset of cell groups by setting sources.use and targets.use.
cellchat <- aggregateNet(cellchat)

# We can also visualize the aggregated cell-cell communication network. For example, showing the number of interactions or the total interaction strength (weights) between any two cell groups using circle plot.
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

# Due to the complicated cell-cell communication network, we can examine the signaling sent from each cell group. 
# Here we also control the parameter edge.weight.max so that we can compare edge weights between different networks.
mat <- H1@net$weight
par(mfrow = c(3,5), xpd=TRUE)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}

pathways.show <- c("LAMININ") 

# Chord diagram
par(mfrow=c(1,1))
pdf(file ="allinteractions.pdf", width = 40, height =40)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()

# Chord diagram with grouping
group.cellType <- c(rep("EMT", 2), rep("stellate", 2), rep("prolif", 2), rep("panc", 5), rep("stem",1), rep("endo",1), rep("toxic",1)) # grouping cell clusters into fibroblast, DC and TC cells
names(group.cellType) <- levels(cellchat@idents)
netVisual_chord_cell(cellchat, signaling = pathways.show, group = group.cellType, title.name = paste0(pathways.show, " signaling network"))

# Heatmap
par(mfrow=c(1,1))
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")

# Show all the significant interactions (L-R pairs) based on user's input (defined by `pairLR.use`)
pairLR.use <- extractEnrichedLR(cellchat, signaling = c("MHC-I"))
netVisual_bubble(cellchat, sources.use = c(1:2), targets.use = c(1), pairLR.use = pairLR.use, remove.isolate = TRUE)

# Compute the network centrality scores
NN3cellchat <- netAnalysis_computeCentrality(NN3cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
# Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
netAnalysis_signalingRole_network(NN3, signaling = pathways.show, width = 8, height = 2.5, font.size = 10)

# Visualize the dominant senders (sources) and receivers (targets) in a 2D space
# We also provide another intutive way to visualize the dominant senders (sources) and receivers (targets) in a 2D space using scatter plot.
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
gg1 <- netAnalysis_signalingRole_scatter(cellchat)
# Signaling role analysis on the cell-cell communication networks of interest
gg2 <- netAnalysis_signalingRole_scatter(cellchat, signaling = c("CXCL", "CCL"))
gg1

# Identify signals contributing most to outgoing or incoming signaling of certain cell groups
# Signaling role analysis on the aggregated cell-cell communication network from all signaling pathways
netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", height = 20, font.size = 6)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", height = 20, font.size = 6, cluster.cols=TRUE)
ht1 + ht2
# Signaling role analysis on the cell-cell communication networks of interest
ht <- netAnalysis_signalingRole_heatmap(cellchat, signaling = c("CXCL", "CCL"))

##### SAVE #####
saveRDS(cellchat, file = "obj_cellchat.rds")

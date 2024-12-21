#-----------------------
# load in 'downsampled_query.rds';
#-----------------------

library(CellChat)
library(patchwork)
library(Seurat)

options(stringsAsFactors = FALSE)

visium.brain<- readRDS('downsampled_query.rds')
metadata<- visium.brain@meta.data

## use annotation results from scPred. 'scpred_no_rejection'
Idents(visium.brain)<- 'scpred_no_rejection'
length(unique(visium.brain$scpred_no_rejection))

visium.brain$scpred_no_rejection2<- ifelse(
  visium.brain$scpred_no_rejection %in% c(
    "Oligo", "VLMC" , "Astro", "SMC-Peri" ,
    "Micro-PVM", "Endo", 'rand'), visium.brain$scpred_no_rejection, 
 'Neuron'
   )

Idents(visium.brain)<- 'scpred_no_rejection2'
length(unique(visium.brain$scpred_no_rejection2))


  
color.use <- scPalette(nlevels(visium.brain)); names(color.use) <- levels(visium.brain)

data.input = Seurat::GetAssayData(visium.brain, slot = "data", assay = "Spatial.008um") # normalized data matrix

meta = data.frame(labels = Seurat::Idents(visium.brain), samples = "sample1", row.names = names(Seurat::Idents(visium.brain))) # manually create a dataframe consisting of the cell labels

meta$samples <- factor(meta$samples)
unique(meta$labels) # check the cell labels
unique(meta$samples) # check the sample labels

spatial.locs = Seurat::GetTissueCoordinates(visium.brain, scale = NULL, cols = c("imagerow", "imagecol")) 

Seurat::SpatialDimPlot(visium.brain, label = T, label.size = 3, cols = color.use)



scalefactors<- jsonlite::fromJSON(txt = file.path("/Users/yining.pan/Desktop/VisiumHD_pipeline/binned_outputs/square_008um/spatial", 'scalefactors_json.json'))
spot.size = 65 # the theoretical spot size (um) in 10X Visium
conversion.factor = spot.size/scalefactors$spot_diameter_fullres
spatial.factors = data.frame(ratio = conversion.factor, tol = spot.size/2)

spatial.locs<- spatial.locs[, -3]

rm(visium.brain); gc();

d.spatial <- computeCellDistance(coordinates = spatial.locs, ratio = spatial.factors$ratio, tol = spatial.factors$tol)
min(d.spatial[d.spatial!=0])



#----------------------
# preprocessing 
#-----------------------
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, spatial.factors = spatial.factors)

CellChatDB <- CellChatDB.mouse # use CellChatDB.human if running on human data
showDatabaseCategory(CellChatDB)
CellChatDB.use <- CellChatDB
cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat, do.fast=F)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = F)


# Part ii: 
## Compute the communication probability and infer cellular communication network
cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1,
                              distance.use = TRUE, interaction.range = 250, scale.distance = 0.01,
                              contact.dependent = TRUE, contact.range = 100)

##Infer the cell-cell communication at a signaling pathway level
cellchat <- computeCommunProbPathway(cellchat)

##Calculate the aggregated cell-cell communication network
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)

netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")

netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")

# part III: Visualization of cell-cell communication network
df.net<- subsetCommunication(cellchat)
pathways.show <- c("WNT") 

# Circle plot
par(mfrow=c(1,1), xpd = TRUE) # `xpd = TRUE` should be added to show the title
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

# Spatial plot
par(mfrow=c(1,1))
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "spatial", edge.width.max = 2, vertex.size.max = 1, alpha.image = 0.2, vertex.label.cex = 3.5)

# Take an input of a few genes
spatialFeaturePlot(cellchat, features = c('Apoe'), point.size = 0.8, color.heatmap = "Reds", direction = 1)


saveRDS(cellchat, file = "cellchat_visiumHD_mouse.rds")
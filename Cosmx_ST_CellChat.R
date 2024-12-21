package.version('CellChat')
##[1] "2.1.2"

#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/CellChat_analysis_of_multiple_spatial_transcriptomics_datasets.html
#https://htmlpreview.github.io/?https://github.com/jinworks/CellChat/blob/master/tutorial/FAQ_on_applying_CellChat_to_spatial_transcriptomics_data.html#cosmx
library(CellChat)
library(patchwork)
library(Seurat)

all.refine.anno<- readRDS('all.hippo.refine.anno.rds')
e3.merge.anno<- subset(all.refine.anno, subset=APOE=='E3')

#-----------------------E3 CRE- -------------------------------------------------------------#
e3.cre_neg_361_1<- subset(e3.merge.anno, subset=sampleID2=='361-1')
e3.cre_neg_361_2<- subset(e3.merge.anno, subset=sampleID2=='361-2')
e3.cre_neg_362_1<- subset(e3.merge.anno, subset=sampleID2=='362-1')
e3.cre_neg_362_2<- subset(e3.merge.anno, subset=sampleID2=='362-2') 

Idents(e3.cre_neg_361_1)<- 'cell_type_me2'
Idents(e3.cre_neg_361_2)<- 'cell_type_me2'
Idents(e3.cre_neg_362_1)<- 'cell_type_me2'
Idents(e3.cre_neg_362_2)<- 'cell_type_me2'

# Prepare input data for CelChat analysis
data.input1 = Seurat::GetAssayData(e3.cre_neg_361_1, slot = "data", assay = "SCT") # normalized data matrix
data.input2 = Seurat::GetAssayData(e3.cre_neg_361_2, slot = "data", assay = "SCT") # normalized data matrix
data.input3 = Seurat::GetAssayData(e3.cre_neg_362_1, slot = "data", assay = "SCT") # normalized data matrix
data.input4 = Seurat::GetAssayData(e3.cre_neg_362_2, slot = "data", assay = "SCT") # normalized data matrix

genes.common <- Reduce(intersect, list(rownames(data.input1), rownames(data.input2), rownames(data.input3),rownames(data.input4) ))

colnames(data.input1) <- paste0("A1_", colnames(data.input1))
colnames(data.input2) <- paste0("A2_", colnames(data.input2))
colnames(data.input3) <- paste0("A3_", colnames(data.input3))
colnames(data.input4) <- paste0("A4_", colnames(data.input4))
data.input <- cbind(data.input1[genes.common, ], data.input2[genes.common, ], 
                    data.input3[genes.common, ], data.input4[genes.common, ])


# define the meta data
# a column named `samples` should be provided for spatial transcriptomics analysis, which is useful for analyzing cell-cell communication by aggregating multiple samples/replicates. Of note, for comparison analysis across different conditions, users still need to create a CellChat object seperately for each condition.  
meta1 = data.frame(labels = Idents(e3.cre_neg_361_1), samples = "A1") # manually create a dataframe consisting of the cell labels
meta2 = data.frame(labels = Idents(e3.cre_neg_361_2), samples = "A2") 
meta3 = data.frame(labels = Idents(e3.cre_neg_362_1), samples = "A3")
meta4 = data.frame(labels = Idents(e3.cre_neg_362_2), samples = "A4")

meta <- rbind(meta1, meta2, meta3, meta4)
rownames(meta) <- colnames(data.input)
# a factor level should be defined for the `meta$labels` and `meta$samples`
meta$labels <- factor(meta$labels, levels = levels(Idents(e3.cre_neg_361_1)))
meta$samples <- factor(meta$samples, levels = c("A1", "A2", "A3", 'A4'))
unique(meta$labels) # check the cell labels

unique(meta$samples) # check the sample labels

spatial.locs1 = Seurat::GetTissueCoordinates(e3.cre_neg_361_1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs2 = Seurat::GetTissueCoordinates(e3.cre_neg_361_2, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs3 = Seurat::GetTissueCoordinates(e3.cre_neg_362_1, scale = NULL, cols = c("imagerow", "imagecol")) 
spatial.locs4 = Seurat::GetTissueCoordinates(e3.cre_neg_362_2, scale = NULL, cols = c("imagerow", "imagecol")) 

spatial.locs1<- spatial.locs1 %>% dplyr::select(1, 2)
spatial.locs2<- spatial.locs2 %>% dplyr::select(1, 2)
spatial.locs3<- spatial.locs3 %>% dplyr::select(1, 2)
spatial.locs4<- spatial.locs4 %>% dplyr::select(1, 2)

spatial.locs <- rbind(spatial.locs1, spatial.locs2,spatial.locs3, spatial.locs4 )
rownames(spatial.locs) <- colnames(data.input)


# Scale factors of spatial coordinates
conversion.factor = 0.18
d1 = computeCellDistance(spatial.locs1)
spot.size = min(d)*conversion.factor # converting the distance in Pixels to Micrometers
spatial.factors1 = data.frame(ratio = conversion.factor, tol = spot.size/2)

conversion.factor = 0.18
d2 = computeCellDistance(spatial.locs2)
spot.size = min(d)*conversion.factor # converting the distance in Pixels to Micrometers
spatial.factors2 = data.frame(ratio = conversion.factor, tol = spot.size/2)

conversion.factor = 0.18
d3 = computeCellDistance(spatial.locs3)
spot.size = min(d)*conversion.factor # converting the distance in Pixels to Micrometers
spatial.factors3 = data.frame(ratio = conversion.factor, tol = spot.size/2)

conversion.factor = 0.18
d4 = computeCellDistance(spatial.locs4)
spot.size = min(d)*conversion.factor # converting the distance in Pixels to Micrometers
spatial.factors4 = data.frame(ratio = conversion.factor, tol = spot.size/2)

spatial.factors <- rbind(spatial.factors1, spatial.factors2, spatial.factors3, spatial.factors4)
rownames(spatial.factors) <- c("A1", "A2", "A3", "A4")

cellchat <- createCellChat(object = data.input, meta = meta, group.by = "labels",
                           datatype = "spatial", coordinates = spatial.locs, 
                           spatial.factors = spatial.factors)
cellchat

CellChatDB <- CellChatDB.mouse
cellchat@DB <- CellChatDB


# subset the expression data of signaling genes for saving computation cost
cellchat <- subsetData(cellchat) # This step is necessary even if using the whole database
future::plan("multisession", workers = 4) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)
#> The number of highly variable ligand-receptor pairs used for signaling inference is 422

contact.range<- mean(min(d1), min(d2), min(d3), min(d4))

cellchat <- computeCommunProb(cellchat, type = "truncatedMean", trim = 0.1, 
                              distance.use = FALSE, interaction.range = 250, scale.distance = NULL,
                              contact.dependent = TRUE, contact.range = contact.range) 
## For spatial transcriptomics in a single-cell resolution, contact.range is approximately equal to the estimated cell diameter (i.e., the cell center-to-center distance), which means that contact-dependent or juxtacrine signaling can only happens when the two cells are contact to each other.  

cellchat <- filterCommunication(cellchat, min.cells = 10)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)


#-----------------
# inferred info |
#-----------------
df.net<- subsetCommunication(cellchat)

unique.LRs<- cellchat@net$LRs

unique.pathways<- cellchat@netP$pathways

write.csv(df.net, 'E3.cre.neg_cellchat_network.csv')

#-----------------
# visualization |
#-----------------

## number of interaction and interaction strength 
groupSize <- as.numeric(table(cellchat@idents))
par(mfrow = c(1,2), xpd=TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = T, label.edge= F, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = T, label.edge= F, title.name = "Interaction weights/strength")
netVisual_heatmap(cellchat, measure = "count", color.heatmap = "Blues")


## certain signaling pathway circle and spatial plot -----------
pathways.show <- c("APP") 

#' Circle plot
par(mfrow=c(1,1), xpd=TRUE)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")

#' Spatial plot
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP") # the slot 'netP' means the inferred intercellular communication network of signaling pathways
par(mfrow=c(1,4))

p1<- netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A1", 
                         layout = "spatial", edge.width.max = 2, alpha.image = 0.2, 
                         vertex.weight = "incoming", vertex.size.max = 6, 
                         vertex.label.cex = 0) & theme(legend.position = "none")
p2<- netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A2", 
                         layout = "spatial", edge.width.max = 2, alpha.image = 0.2, 
                         vertex.weight = "incoming", vertex.size.max = 6, 
                         vertex.label.cex = 0)& theme(legend.position = "none")
p3<- netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A3", 
                         layout = "spatial", edge.width.max = 2, alpha.image = 0.2, 
                         vertex.weight = "incoming", vertex.size.max = 6, 
                         vertex.label.cex = 0)& theme(legend.position = "none")
p4<- netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A4", 
                         layout = "spatial", edge.width.max = 2, alpha.image = 0.2, 
                         vertex.weight = "incoming", vertex.size.max = 6, 
                         vertex.label.cex = 0)
library(patchwork)
pdf('E3_neg_four_slides_app_signaling_on_spatial_circle_by_incoming.pdf', width = 20, height = 10)
(p1 | p2 | p3 | p4)
graphics.off()


p1<- netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A1", 
                         layout = "spatial", edge.width.max = 2, alpha.image = 0.2, 
                         vertex.weight = "outgoing", vertex.size.max = 6, 
                         vertex.label.cex = 0) & theme(legend.position = "none")
p2<- netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A2", 
                         layout = "spatial", edge.width.max = 2, alpha.image = 0.2, 
                         vertex.weight = "outgoing", vertex.size.max = 6, 
                         vertex.label.cex = 0)& theme(legend.position = "none")
p3<- netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A3", 
                         layout = "spatial", edge.width.max = 2, alpha.image = 0.2, 
                         vertex.weight = "outgoing", vertex.size.max = 6, 
                         vertex.label.cex = 0)& theme(legend.position = "none")
p4<- netVisual_aggregate(cellchat, signaling = pathways.show, sample.use = "A4", 
                         layout = "spatial", edge.width.max = 2, alpha.image = 0.2, 
                         vertex.weight = "outgoing", vertex.size.max = 6, 
                         vertex.label.cex = 0)
library(patchwork)
pdf('E3_neg_four_slides_app_signaling_on_spatial_circle_by_outgoing.pdf', width = 20, height = 10)
(p1 | p2 | p3 | p4)
graphics.off()


## Visualize the computed centrality scores using heatmap, allowing ready identification of major signaling roles of cell groups
par(mfrow=c(1,1))
netAnalysis_signalingRole_network(cellchat, signaling = pathways.show, 
                                  width = 8, height = 2.5, font.size = 10)


## Compute the contribution of each LR to the overall signaling pathway ---------
netAnalysis_contribution(cellchat, signaling = pathways.show)


## Take an input of a few genes -----------
p1<- spatialFeaturePlot(cellchat, features = c("App","Trem2", "Cd74", "Tyrobp"), 
                        sample.use = "A1", point.size = 0.8, 
                        color.heatmap = "Reds", direction = 1)
p2<- spatialFeaturePlot(cellchat, features = c("App","Trem2", "Cd74", "Tyrobp"), 
                        sample.use = "A2", point.size = 0.8, 
                        color.heatmap = "Reds", direction = 1)
p3<- spatialFeaturePlot(cellchat, features = c("App","Trem2", "Cd74", "Tyrobp"), 
                        sample.use = "A3", point.size = 0.8, 
                        color.heatmap = "Reds", direction = 1)
p4<- spatialFeaturePlot(cellchat, features = c("App","Trem2", "Cd74", "Tyrobp"), 
                        sample.use = "A4", point.size = 0.8, 
                        color.heatmap = "Reds", direction = 1)
library(patchwork)
pdf('E3_neg_four_slides_app_signaling_genes_expression.pdf', width = 20, height = 10)
(p1 / p2 / p3 / p4)
graphics.off()

## Take an input of a ligand-receptor pair ------------

p1<- spatialFeaturePlot(cellchat, pairLR.use = "APP_TREM2_TYROBP", sample.use = "A1", 
                        point.size = 1.5, do.binary = TRUE, cutoff = 0.05, 
                        enriched.only = F, color.heatmap = "Reds", direction = 1)
p2<- spatialFeaturePlot(cellchat, pairLR.use = "APP_TREM2_TYROBP", sample.use = "A2", 
                        point.size = 1.5, do.binary = TRUE, cutoff = 0.05, 
                        enriched.only = F, color.heatmap = "Reds", direction = 1)
p3<- spatialFeaturePlot(cellchat, pairLR.use = "APP_TREM2_TYROBP", sample.use = "A3", 
                        point.size = 1.5, do.binary = TRUE, cutoff = 0.05, 
                        enriched.only = F, color.heatmap = "Reds", direction = 1)
p4<- spatialFeaturePlot(cellchat, pairLR.use = "APP_TREM2_TYROBP", sample.use = "A4", 
                        point.size = 1.5, do.binary = TRUE, cutoff = 0.05, 
                        enriched.only = F, color.heatmap = "Reds", direction = 1)

library(patchwork)
pdf('E3_neg_four_slides_app_signaling_genes_colocalize_on_spatial.pdf', width = 20, height = 10)
(p1 | p2 | p3 | p4)
graphics.off()

## save for later use

saveRDS(cellchat_e3_neg, file = "cellchat_MOUSE_E3_CRE_neg_ReplicatesA1A2A3A4.rds")

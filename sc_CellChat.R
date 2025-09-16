
library(CellChat)

integrated.anno.sub<- readRDS('integrated.anno.sub.rds')


#split the integrated.anno into p56 and GFP
integrated.anno.sub.p56<- subset(integrated.anno.sub, subset=genotype == 'AAV-p56')
integrated.anno.sub.GFP<- subset(integrated.anno.sub, subset=genotype == 'AAV-GFP')

#MAKE cellchat object for p56 and GFP
cellChat.p56 <- createCellChat(object = integrated.anno.sub.p56, group.by = "ident", assay = "RNA")
cellChat.GFP <- createCellChat(object = integrated.anno.sub.GFP, group.by = "ident", assay = "RNA")


#analyze each one
levels(cellChat.p56@idents)
levels(cellChat.GFP@idents)


CellChatDB <- CellChatDB.mouse
CellChatDB.use <- CellChatDB

cellChat.p56@DB <- CellChatDB.use
cellChat.GFP@DB <- CellChatDB.use


cellChat.p56 <- subsetData(cellChat.p56)
future::plan("multisession", workers = 4) # do parallel
cellChat.p56 <- identifyOverExpressedGenes(cellChat.p56)
cellChat.p56 <- identifyOverExpressedInteractions(cellChat.p56)


cellChat.GFP <- subsetData(cellChat.GFP)
future::plan("multisession", workers = 4) # do parallel
cellChat.GFP <- identifyOverExpressedGenes(cellChat.GFP)
cellChat.GFP <- identifyOverExpressedInteractions(cellChat.GFP)


cellChat.p56 <- computeCommunProb(cellChat.p56, type = "triMean")
cellChat.GFP <- computeCommunProb(cellChat.GFP, type = "triMean")

cellChat.p56 <- computeCommunProbPathway(cellChat.p56)
cellChat.GFP <- computeCommunProbPathway(cellChat.GFP)

cellChat.p56 <- aggregateNet(cellChat.p56)
cellChat.GFP <- aggregateNet(cellChat.GFP)

# Compute the network centrality scores
cellChat.p56 <- netAnalysis_computeCentrality(cellChat.p56, slot.name = "netP")
cellChat.GFP <- netAnalysis_computeCentrality(cellChat.GFP, slot.name = "netP")


#filter out the cell-cell communication with only few cells in certain cell groups.
cellChat.p56 <- filterCommunication(cellChat.p56, min.cells = 10)
cellChat.GFP <- filterCommunication(cellChat.GFP, min.cells = 10)


#merge the two 
object.list <- list(GFP = cellChat.GFP, p56 = cellChat.p56)
cellchat <- mergeCellChat(object.list, add.names = names(object.list))
cellchat

##save for later use
save(object.list, file = "cellchat_object.list_kj_mouse_GFP_p56_sub.RData")
save(cellchat, file = "cellchat_merged_kj_mouse_GFP_p56_sub.RData")



#---------------------
# Visualization
#---------------------

##define a function to remove cell labels on circle plots showing diff #interactions and strength
netVisual_diffInteraction_noLabel <- function (object, comparison = c(1, 2), measure = c("count", 
                                                                                         "weight", "count.merged", "weight.merged"), color.use = NULL, 
                                               color.edge = c("#b2182b", "#2166ac"), title.name = NULL, 
                                               sources.use = NULL, targets.use = NULL, remove.isolate = FALSE, 
                                               top = 1, weight.scale = FALSE, vertex.weight = 20, vertex.weight.max = NULL, 
                                               vertex.size.max = 15, 
                                               edge.weight.max = NULL, edge.width.max = 8, alpha.edge = 0.6, 
                                               label.edge = FALSE, 
                                               edge.curved = 0.2, shape = "circle", layout = in_circle(), 
                                               margin = 0.2, arrow.width = 1, arrow.size = 0.2) 
{
  options(warn = -1)
  measure <- match.arg(measure)
  obj1 <- object@net[[comparison[1]]][[measure]]
  obj2 <- object@net[[comparison[2]]][[measure]]
  net.diff <- obj2 - obj1
  if (measure %in% c("count", "count.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential number of interactions"
    }
  }
  else if (measure %in% c("weight", "weight.merged")) {
    if (is.null(title.name)) {
      title.name = "Differential interaction strength"
    }
  }
  net <- net.diff
  if ((!is.null(sources.use)) | (!is.null(targets.use))) {
    df.net <- reshape2::melt(net, value.name = "value")
    colnames(df.net)[1:2] <- c("source", "target")
    if (!is.null(sources.use)) {
      if (is.numeric(sources.use)) {
        sources.use <- rownames(net.diff)[sources.use]
      }
      df.net <- subset(df.net, source %in% sources.use)
    }
    if (!is.null(targets.use)) {
      if (is.numeric(targets.use)) {
        targets.use <- rownames(net.diff)[targets.use]
      }
      df.net <- subset(df.net, target %in% targets.use)
    }
    cells.level <- rownames(net.diff)
    df.net$source <- factor(df.net$source, levels = cells.level)
    df.net$target <- factor(df.net$target, levels = cells.level)
    df.net$value[is.na(df.net$value)] <- 0
    net <- tapply(df.net[["value"]], list(df.net[["source"]], 
                                          df.net[["target"]]), sum)
    net[is.na(net)] <- 0
  }
  if (remove.isolate) {
    idx1 <- which(Matrix::rowSums(net) == 0)
    idx2 <- which(Matrix::colSums(net) == 0)
    idx <- intersect(idx1, idx2)
    net <- net[-idx, ]
    net <- net[, -idx]
  }
  net[abs(net) < stats::quantile(abs(net), probs = 1 - top, 
                                 na.rm = T)] <- 0
  g <- graph_from_adjacency_matrix(net, mode = "directed", 
                                   weighted = T)
  edge.start <- igraph::ends(g, es = igraph::E(g), names = FALSE)
  coords <- layout_(g, layout)
  if (nrow(coords) != 1) {
    coords_scale = scale(coords)
  }
  else {
    coords_scale <- coords
  }
  if (is.null(color.use)) {
    color.use = scPalette(length(igraph::V(g)))
  }
  if (is.null(vertex.weight.max)) {
    vertex.weight.max <- max(vertex.weight)
  }
  vertex.weight <- vertex.weight/vertex.weight.max * vertex.size.max + 
    5
  loop.angle <- ifelse(coords_scale[igraph::V(g), 1] > 0, -atan(coords_scale[igraph::V(g), 
                                                                             2]/coords_scale[igraph::V(g), 1]), pi - atan(coords_scale[igraph::V(g), 
                                                                                                                                       2]/coords_scale[igraph::V(g), 1]))
  igraph::V(g)$size <- vertex.weight
  igraph::V(g)$color <- color.use[igraph::V(g)]
  igraph::V(g)$frame.color <- color.use[igraph::V(g)]
  if (label.edge) {
    igraph::E(g)$label <- igraph::E(g)$weight
    igraph::E(g)$label <- round(igraph::E(g)$label, digits = 1)
  }
  igraph::E(g)$arrow.width <- arrow.width
  igraph::E(g)$arrow.size <- arrow.size
  igraph::E(g)$color <- ifelse(igraph::E(g)$weight > 0, color.edge[1], 
                               color.edge[2])
  igraph::E(g)$color <- grDevices::adjustcolor(igraph::E(g)$color, 
                                               alpha.edge)
  igraph::E(g)$weight <- abs(igraph::E(g)$weight)
  if (is.null(edge.weight.max)) {
    edge.weight.max <- max(igraph::E(g)$weight)
  }
  if (weight.scale == TRUE) {
    igraph::E(g)$width <- 0.3 + igraph::E(g)$weight/edge.weight.max * 
      edge.width.max
  }
  else {
    igraph::E(g)$width <- 0.3 + edge.width.max * igraph::E(g)$weight
  }
  if (sum(edge.start[, 2] == edge.start[, 1]) != 0) {
    igraph::E(g)$loop.angle[which(edge.start[, 2] == edge.start[, 
                                                                1])] <- loop.angle[edge.start[which(edge.start[, 2] == edge.start[, 1]), 1]]
  }
  radian.rescale <- function(x, start = 0, direction = 1) {
    c.rotate <- function(x) (x + start)%%(2 * pi) * direction
    c.rotate(scales::rescale(x, c(0, 2 * pi), range(x)))
  }
  plot(g, edge.curved = edge.curved, vertex.shape = shape, 
       layout = coords_scale, margin = margin,
       vertex.label = NA,  # Disable vertex labels
       edge.label = NA      # Disable edge labels
  )
  if (!is.null(title.name)) {
    text(0, 1.5, title.name, cex = 1.1)
  }
  gg <- recordPlot()
  return(gg)
}


#bar charts showing #interactions & #intensities------------------------------------------
gg1 <- compareInteractions(cellchat, show.legend = F, group = c(1,2))
gg2 <- compareInteractions(cellchat, show.legend = F, group = c(1,2), measure = "weight")
gg1 + gg2


#circle plots showing #interactions & #intensities----------------------
par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction_noLabel(cellchat, weight.scale = T)
netVisual_diffInteraction_noLabel(cellchat, weight.scale = T, measure = "weight")

par(mfrow = c(1,2), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")



#CCI players-------------------------------------------------------------------
num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, do.label = F)
}


x_limits <- c(min(sapply(gg, function(g) min(g$data$x))), max(sapply(gg, function(g) max(g$data$x))))
y_limits <- c(min(sapply(gg, function(g) min(g$data$y))), max(sapply(gg, function(g) max(g$data$y))))
gg_modified <- list()

for (i in 1:length(gg)) {
  gg_modified[[i]] <- gg[[i]] + xlim(x_limits) + ylim(y_limits)
}

pdf('major_cci_sub_players.pdf',height=6, width=10)
patchwork::wrap_plots(plots = gg_modified)
graphics.off()


#Information workflow---------------------------------------------------------
pdf('information_workflow_sub.pdf', width = 14,height = 8)
print(rankNet(cellchat, mode = "comparison", measure = "weight", sources.use = NULL, targets.use = NULL, stacked = F, do.stat = TRUE))
graphics.off()


#outgoing & incoming signaling patterns---------------------------------------------------------
library(ComplexHeatmap)
i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "outgoing", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15)
pdf('outgoing_signaling patterns_sub.pdf', width = 14,height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
graphics.off()


i = 1
pathway.union <- union(object.list[[i]]@netP$pathways, object.list[[i+1]]@netP$pathways)
ht1 = netAnalysis_signalingRole_heatmap(object.list[[i]],pattern = "incoming", signaling = pathway.union, title = names(object.list)[i], width = 5, height = 15)
ht2 = netAnalysis_signalingRole_heatmap(object.list[[i+1]], pattern = "incoming", signaling = pathway.union, title = names(object.list)[i+1], width = 5, height = 15)
pdf('incoming_signaling patterns_sub.pdf', width = 14,height = 15)
draw(ht1 + ht2, ht_gap = unit(0.5, "cm"))
graphics.off()

# LR pairs changes, 1. Identify dysfunctional signaling by comparing the communication probabities----------------------------------------------------------------------
#' This method for identifying the upgulated and down-regulated signaling is perfomed by 
#' comparing the communication probability between two datasets for each L-R pair and each pair of cell groups.

###venous--------------------
gg1 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:6), comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in p56", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = 1, targets.use = c(2:6),  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in p56", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object

pdf('ligand-pairs changes_venous_as_source.pdf', width = 14,height = 15)
gg1 + gg2
graphics.off()


gg1 <- netVisual_bubble(cellchat, sources.use = c(2:6), targets.use = 1,  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in p56", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object
gg2 <- netVisual_bubble(cellchat, sources.use = c(2:6), targets.use = 1,  comparison = c(1, 2), max.dataset = 1, title.name = "Decreased signaling in p56", angle.x = 45, remove.isolate = T)
#> Comparing communications on a merged object

pdf('ligand-pairs changes_venous_as_target.pdf', width = 14,height = 15)
gg1 + gg2
graphics.off()

# LR pairs changes, 2. Identify dysfunctional signaling by DEGs--------------------------------------------------------------------------------------------------------------
#' This alternative method is to identify the upgulated and down-regulated signaling ligand-receptor pairs based on the differential expression analysis (DEA). 
#' Specifically, we perform differential expression analysis between two biological conditions (i.e., NL and LS) for each cell group, 
#' and then obtain the upgulated and down-regulated signaling based on the fold change of ligands in the sender cells and receptors in the receiver cells.
#' Of note, users may observe the same LR pairs appearing in both the up-regulated and down-regulated results due to the fact that DEA between conditions is performed for each cell group. 
#' There is a param to allow users to perform DEA between conditions by ignoring cell group information


pos.dataset = "ko"
features.name = paste0(pos.dataset, ".merged")
cellchat <- identifyOverExpressedGenes(cellchat, group.dataset = "datasets",
                                       pos.dataset = pos.dataset, features.name = features.name, 
                                       only.pos = FALSE, thresh.pc = 0.1, thresh.fc = 0.05,thresh.p = 0.05, group.DE.combined = FALSE) 

# map the results of differential expression analysis onto the inferred cell-cell communications to easily manage/subset the ligand-receptor pairs of interest
net <- netMappingDEG(cellchat, features.name = features.name, variable.all = TRUE)
# extract the ligand-receptor pairs with upregulated ligands in ko
net.up <- subsetCommunication(cellchat, net = net, datasets = "ko",ligand.logFC = 0.05, receptor.logFC = NULL)
# extract the ligand-receptor pairs with upregulated ligands and upregulated receptors in wt, i.e.,downregulated in ko
net.down <- subsetCommunication(cellchat, net = net, datasets = "wt",ligand.logFC = -0.05, receptor.logFC = NULL)

gene.up <- extractGeneSubsetFromPair(net.up, cellchat)
gene.down <- extractGeneSubsetFromPair(net.down, cellchat)

# Chord diagram
par(mfrow = c(1,2), xpd=TRUE)
netVisual_chord_gene(object.list[[2]], sources.use = 1, targets.use = c(2:3), 
                     slot.name = 'net', net = net.up, lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Up-regulated signaling in ", names(object.list)[2]),
                     color.use = col.map)
netVisual_chord_gene(object.list[[1]], sources.use = 1, targets.use = c(2:3), slot.name = 'net',
                     net = net.down %>% filter(source == 'ExN1', target %in% c('ExN1', 'ExN2', 'ExN3')) %>% slice_head(n=50),  # here im just subsetting the net.down dataframe as the plot shows too many LRs
                     lab.cex = 0.8, small.gap = 3.5, 
                     title.name = paste0("Down-regulated signaling in ", names(object.list)[2]),
                     color.use = col.map, show.legend = F)

write.csv(net.down, 's192_LR_changes_DEG_net.down.csv')
write.csv(net.up, 's192_LR_changes_DEG_net.up.csv')


                                                                     
##explore LR pairs
df.net.GFP <- subsetCommunication(cellChat.GFP)
df.net.p56 <- subsetCommunication(cellChat.p56)



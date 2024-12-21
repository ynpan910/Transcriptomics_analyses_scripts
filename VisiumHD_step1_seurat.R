library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)


#------------------
# load Visium HD data
#------------------

## download from here https://www.10xgenomics.com/datasets/visium-hd-cytassist-gene-expression-libraries-of-mouse-brain-he

localdir <- "/Users/yining.pan/Desktop/VisiumHD_pipeline/"
list.files(localdir)
object <- Load10X_Spatial(data.dir = localdir, bin.size = c(8, 16))

# Setting default assay changes between 8um and 16um binning
Assays(object)
DefaultAssay(object) <- "Spatial.008um"

vln.plot <- VlnPlot(object, features = "nCount_Spatial.008um", pt.size = 0) + theme(axis.text = element_text(size = 4)) + NoLegend()
count.plot <- SpatialFeaturePlot(object, features = "nCount_Spatial.008um") + theme(legend.position = "right")
vln.plot | count.plot


#------------------
# normalize datasets
#------------------

# normalize both 8um and 16um bins
DefaultAssay(object) <- "Spatial.008um"
object <- NormalizeData(object)

DefaultAssay(object) <- "Spatial.016um"
object <- NormalizeData(object)



#-------------------------
# visualize gene expression
#-------------------------
# switch spatial resolution to 16um from 8um
DefaultAssay(object) <- "Spatial.016um"
p1 <- SpatialFeaturePlot(object, features = "Rorb") + ggtitle("Rorb expression (16um)")

# switch back to 8um
DefaultAssay(object) <- "Spatial.008um"
p2 <- SpatialFeaturePlot(object, features = "Hpca") + ggtitle("Hpca expression (8um)")

p1 | p2


#-------------------------
# Unsupervised clustering
#-------------------------

# note that data is already normalized
DefaultAssay(object) <- "Spatial.008um"
object <- FindVariableFeatures(object)
object <- ScaleData(object)
# we select 50,0000 cells and create a new 'sketch' assay
object <- SketchData(
  object = object,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

# switch analysis to sketched cells
DefaultAssay(object) <- "sketch"

# perform clustering workflow
object <- FindVariableFeatures(object)
object <- ScaleData(object)
object <- RunPCA(object, assay = "sketch", reduction.name = "pca.sketch")
object <- FindNeighbors(object, assay = "sketch", reduction = "pca.sketch", dims = 1:50)
object <- FindClusters(object, cluster.name = "seurat_cluster.sketched", resolution = 3)
object <- RunUMAP(object, reduction = "pca.sketch", reduction.name = "umap.sketch", return.model = T, dims = 1:50)


object <- ProjectData(
  object = object,
  assay = "Spatial.008um",
  full.reduction = "full.pca.sketch",
  sketched.assay = "sketch",
  sketched.reduction = "pca.sketch",
  umap.model = "umap.sketch",
  dims = 1:50,
  refdata = list(seurat_cluster.projected = "seurat_cluster.sketched")
)


#We can visualize the clustering results for the sketched cells, as well as the projected clustering results for the full dataset:
DefaultAssay(object) <- "sketch"
Idents(object) <- "seurat_cluster.sketched"
p1 <- DimPlot(object, reduction = "umap.sketch", label = F) + ggtitle("Sketched clustering (50,000 cells)") + theme(legend.position = "bottom")

# switch to full dataset
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
p2 <- DimPlot(object, reduction = "full.umap.sketch", label = F) + ggtitle("Projected clustering (full dataset)") + theme(legend.position = "bottom")

p1 | p2

# visualize the clusters on spatial location
SpatialDimPlot(object, label = T, repel = T, label.size = 4)

# plot a cluster at a time
Idents(object) <- "seurat_cluster.projected"
cells <- CellsByIdentities(object, idents = c(0, 4, 32, 34, 35))
p <- SpatialDimPlot(object,
                    cells.highlight = cells[setdiff(names(cells), "NA")],
                    cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T
) + NoLegend()
p

## We can also find and visualize the top gene expression markers for each cluster:
  
# Crete downsampled object to make visualization either
DefaultAssay(object) <- "Spatial.008um"
Idents(object) <- "seurat_cluster.projected"
object_subset <- subset(object, cells = Cells(object[["Spatial.008um"]]), downsample = 1000)

saveRDS(object, 'object.rds')

# Order clusters by similarity
DefaultAssay(object_subset) <- "Spatial.008um"
Idents(object_subset) <- "seurat_cluster.projected"
object_subset <- BuildClusterTree(object_subset, assay = "Spatial.008um", reduction = "full.pca.sketch", reorder = T)

markers <- FindAllMarkers(object_subset, assay = "Spatial.008um", only.pos = TRUE)
markers %>%
  group_by(cluster) %>%
  dplyr::filter(avg_log2FC > 1) %>%
  slice_head(n = 5) %>%
  ungroup() -> top5

object_subset <- ScaleData(object_subset, assay = "Spatial.008um", features = top5$gene)
p <- DoHeatmap(object_subset, assay = "Spatial.008um", features = top5$gene, size = 2.5) + theme(axis.text = element_text(size = 5.5)) + NoLegend()
p



#----------------------------------------------
# Identifying spatially-defined tissue domains
#----------------------------------------------
if (!requireNamespace("Banksy", quietly = TRUE)) {
  remotes::install_github("prabhakarlab/Banksy@devel")
}
library(SeuratWrappers)
library(Banksy)
library(Seurat)

object<- readRDS('object.rds')
object <- RunBanksy(object,
                    lambda = 0.8, verbose = TRUE,
                    assay = "Spatial.008um", slot = "data", features = "variable",
                    k_geom = 50
)



DefaultAssay(object) <- "BANKSY"
object <- RunPCA(object, assay = "BANKSY", reduction.name = "pca.banksy", features = rownames(object), npcs = 30)
object <- FindNeighbors(object, reduction = "pca.banksy", dims = 1:30)
object <- FindClusters(object, cluster.name = "banksy_cluster", resolution = 0.5)

Idents(object) <- "banksy_cluster"
p <- SpatialDimPlot(object, group.by = "banksy_cluster", label = T, repel = T, label.size = 4, , images="slice1.008um")
p

# highlight the spatial location of each tissue domain individually:
banksy_cells <- CellsByIdentities(object)
p <- SpatialDimPlot(object, cells.highlight = banksy_cells[setdiff(names(banksy_cells), "NA")], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T) + NoLegend()
p


saveRDS(object, 'object1.rds')
#--------------------------------
# Subset out anatomical regions |
#--------------------------------
library(Seurat)
library(ggplot2)
library(patchwork)
library(dplyr)
library(sp)

cortex.coordinates <- as.data.frame(read.csv("cortex-hippocampus_coordinates.csv"))
cortex <- CreateSegmentation(cortex.coordinates)

object[["cortex"]] <- Overlay(object[["slice1.008um"]], cortex)
cortex <- subset(object, cells = Cells(object[["cortex"]]))

rm(object)
gc()

#-----------------------------------------
# Integration with scRNA-seq data (deconvolution) |
#-----------------------------------------


# sketch the cortical subset of the Visium HD dataset
DefaultAssay(cortex) <- "Spatial.008um"
cortex <- FindVariableFeatures(cortex)
cortex <- SketchData(
  object = cortex,
  ncells = 50000,
  method = "LeverageScore",
  sketched.assay = "sketch"
)

DefaultAssay(cortex) <- "sketch"
cortex <- ScaleData(cortex)
cortex <- RunPCA(cortex, assay = "sketch", reduction.name = "pca.cortex.sketch", verbose = T)
cortex <- FindNeighbors(cortex, reduction = "pca.cortex.sketch", dims = 1:50)
cortex <- RunUMAP(cortex, reduction = "pca.cortex.sketch", reduction.name = "umap.cortex.sketch", return.model = T, dims = 1:50, verbose = T)

# load in the reference scRNA-seq dataset
ref <- readRDS("allen_scRNAseq_ref.Rds")

Idents(ref) <- "subclass_label"
counts <- ref[["RNA"]]$counts
cluster <- as.factor(ref$subclass_label)
nUMI <- ref$nCount_RNA
levels(cluster) <- gsub("/", "-", levels(cluster))
cluster <- droplevels(cluster)

# create the RCTD reference object
library(spacexr)
reference <- Reference(counts, cluster, nUMI)

counts_hd <- cortex[["sketch"]]$counts
cortex_cells_hd <- colnames(cortex[["sketch"]])
coords <- GetTissueCoordinates(cortex)[cortex_cells_hd, 1:2]

# create the RCTD query object
query <- SpatialRNA(coords, counts_hd, colSums(counts_hd))

rm(counts_hd)
rm(counts)
rm(ref)
gc()

# run RCTD
RCTD <- create.RCTD(query, reference, max_cores = 28)
RCTD <- run.RCTD(RCTD, doublet_mode = "doublet")
saveRDS(RCTD, 'RCTD.rds')

# add results back to Seurat object
cortex <- AddMetaData(cortex, metadata = RCTD@results$results_df)

saveRDS(reference,'reference.rds'); saveRDS(cortex, 'cortex.rds')
rm(reference); gc()
# project RCTD labels from sketched cortical cells to all cortical cells
cortex$first_type <- as.character(cortex$first_type)
cortex$first_type[is.na(cortex$first_type)] <- "Unknown"
cortex <- ProjectData(
  object = cortex,
  assay = "Spatial.008um",
  full.reduction = "pca.cortex",
  sketched.assay = "sketch",
  sketched.reduction = "pca.cortex.sketch",
  umap.model = "umap.cortex.sketch",
  dims = 1:50,
  refdata = list(full_first_type = "first_type")
)

object<- readRDS('object1.rds')
DefaultAssay(object) <- "Spatial.008um"

# we only ran RCTD on the cortical cells
# set labels to all other cells as "Unknown"
object[[]][, "full_first_type"] <- "Unknown"
object$full_first_type[Cells(cortex)] <- cortex$full_first_type[Cells(cortex)]


Idents(object) <- "full_first_type"

# now we can spatially map the location of any scRNA-seq cell type
# start with Layered (starts with L), excitatory neurons in the cortex
cells <- CellsByIdentities(object)
excitatory_names <- sort(grep("^L.* CTX", names(cells), value = TRUE))
p <- SpatialDimPlot(object, cells.highlight = cells[excitatory_names], cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T, ncol = 4)
p


# We can now look for associations between the scRNA-seq labels of individual bins, and their tissue domain identity (as assigned by BANKSY). By asking which domains the excitatory neuron cells fall in, we can rename the BANKSY clusters as neuronal layers:

plot_cell_types <- function(data, label) {
  p <- ggplot(data, aes(x = get(label), y = n, fill = full_first_type)) +
    geom_bar(stat = "identity", position = "stack") +
    geom_text(aes(label = ifelse(n >= min_count_to_show_label, full_first_type, "")), position = position_stack(vjust = 0.5), size = 2) +
    xlab(label) +
    ylab("# of Spots") +
    ggtitle(paste0("Distribution of Cell Types across ", label)) +
    theme_minimal()
}

cell_type_banksy_counts <- object[[]] %>%
  dplyr::filter(full_first_type %in% excitatory_names) %>%
  dplyr::count(full_first_type, banksy_cluster)

min_count_to_show_label <- 20

p <- plot_cell_types(cell_type_banksy_counts, "banksy_cluster")
p

#Based on this plot, we can now assign cells (even if they are not excitatory neurons) to individual neuronal layers.

Idents(object) <- "banksy_cluster"
object$layer_id <- "Unknown"
object$layer_id[WhichCells(object, idents = c(6))] <- "Layer 2/3"
object$layer_id[WhichCells(object, idents = c(11))] <- "Layer 4/5"
object$layer_id[WhichCells(object, idents = c(2))] <- "Layer 6"

# set ID to RCTD label
Idents(object) <- "full_first_type"

# Visualize distribution of 4 interneuron subtypes
inhibitory_names <- c("Sst", "Pvalb", "Vip", "Lamp5")
cell_ids <- CellsByIdentities(object, idents = inhibitory_names)
p <- SpatialDimPlot(object, cells.highlight = cell_ids, cols.highlight = c("#FFFF00", "grey50"), facet.highlight = T, combine = T, ncol = 4)
p

# create barplot to show proportions of cell types of interest
layer_table <- table(object$full_first_type, object$layer_id)[inhibitory_names, 1:4]

neuron_props <- reshape2::melt(prop.table(layer_table), margin = 1)
ggplot(neuron_props, aes(x = Var1, y = value, fill = Var2)) +
  geom_bar(stat = "identity", position = "fill") +
  labs(x = "Cell type", y = "Proportion", fill = "Layer") +
  theme_classic()

saveRDS(object, 'object.final.rds')

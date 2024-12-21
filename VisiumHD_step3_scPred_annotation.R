library("scPred")
library("Seurat")
library("magrittr")

#----------------------
# load in original data
#----------------------
load('Seurat.ss.rda')
reference <- ss.seurat
rm(ss.seurat); gc()
query <- readRDS('object.final.rds')

#----------------------
# downsize both reference and query
#----------------------
#ref
metadata_ref<- reference@meta.data

subclass_labels <- unique(reference$subclass_label)
selected_cells <- c()

for (label in subclass_labels) {
  # Get the cells for the current subclass label
  subclass_cells <- WhichCells(reference, expression = subclass_label == label)
  
  # If there are more than 100 cells, sample 100, otherwise keep all
  if (length(subclass_cells) > 100) {
    subclass_cells <- sample(subclass_cells, 100)
  }
  
  # Add the selected cells for the current subclass to the vector
  selected_cells <- c(selected_cells, subclass_cells)
}

downsampled_reference <- subset(reference, cells = selected_cells)

rm(reference); rm(metadata_ref); gc()

# spatial seurat (query)
metadata_query<- query@meta.data

subclass_labels <- unique(query$category_SingleCellNet)
selected_cells <- c()

for (label in subclass_labels) {
  # Get the cells for the current subclass label
  subclass_cells <- WhichCells(query, expression = category_SingleCellNet == label)
  
  # If there are more than 100 cells, sample 100, otherwise keep all
  if (length(subclass_cells) > 100) {
    subclass_cells <- sample(subclass_cells, 100)
  }
  
  # Add the selected cells for the current subclass to the vector
  selected_cells <- c(selected_cells, subclass_cells)
}

downsampled_query <- subset(query, cells = selected_cells)

rm(query); rm(metadata_query); gc()


reference<- downsampled_reference
query<- downsampled_query

#---------------
# start scPred
#---------------
reference <- reference %>% 
  NormalizeData() %>% 
  FindVariableFeatures() %>% 
  ScaleData() %>% 
  RunPCA() %>% 
  RunUMAP(dims = 1:30)

DimPlot(reference, group.by = "subclass_label", label = TRUE, repel = TRUE)

reference <- getFeatureSpace(reference, "subclass_label")

reference <- trainModel(reference)

get_probabilities(reference) %>% head()

get_scpred(reference)

query <- NormalizeData(query)


## need to change the assay name of Spatial0.08 to RNA, just to meet the scPred's coding
query1<- query
query[['RNA']] = query[['Spatial.008um']]


query <- scPredict(query, reference)

DimPlot(query, group.by = "scpred_prediction", reduction = "scpred")

metadata_query<- query@meta.data


saveRDS(query, 'downsampled_query.rds')
saveRDS(reference, 'downsampled_reference.rds')

## the annotation results are saved in "scpred_prediction"       &       "scpred_no_rejection"  in metadata
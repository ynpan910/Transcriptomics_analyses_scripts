#--------------this is doing infectious samples--------------#

#'tutorial https://cole-trapnell-lab.github.io/monocle-release/docs/#estimate-size-factors-and-dispersions-required
#'https://bookdown.org/ytliu13207/SingleCellMultiOmicsDataAnalysis/monocle2.html#monocle2-process
#'https://cole-trapnell-lab.github.io/monocle-release/tutorial_data/Olsson_dataset_analysis_final.html



library(monocle) #"2.32.0"
library(tidyverse)
library(Seurat)


integrated_inf<- readRDS('integrated_inf.rds')

DefaultAssay(integrated_inf)<- 'RNA'
integrated_inf[["RNA"]] <- JoinLayers(integrated_inf[["RNA"]])
cds<- as.CellDataSet(integrated_inf)
cds@phenoData@data$stage <- factor(cds@phenoData@data$stage, levels = c('Early', 'Medium', 'Late'))
metadata<- cds@phenoData@data


# Estimate size factors and dispersions 
cds <- estimateSizeFactors(cds)
cds <- estimateDispersions(cds)

# Filtering low-quality cells
cds <- detectGenes(cds , min_expr = 0.1)
print(head(fData(cds )))
expressed_genes <- row.names(subset(fData(cds),
                                    num_cells_expressed >= 3))

# unsupervised clustering
disp_table <- dispersionTable(cds)
unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
cds <- setOrderingFilter(cds, unsup_clustering_genes$gene_id)

#plot_ordering_genes(cds)
#plot_pc_variance_explained(cds, return_all = F)

cds <- reduceDimension(cds, max_components = 2, num_dim = 6,
                        reduction_method = 'tSNE', verbose = F) #,
                       #residualModelFormulaStr = "~num_genes_expressed")
cds  <- clusterCells(cds, num_clusters = 2 )

plot_cell_clusters(cds, 1, 2, color_by = 'sampleID')

# constructing single cell trajectories

## Trajectory step 1: choosing genes that define progress
diff_test_res <- differentialGeneTest(cds[expressed_genes,],
                                      fullModelFormulaStr = '~stage')
ordering_genes <- row.names (subset(diff_test_res, qval < 0.01))

cds <- setOrderingFilter(cds, ordering_genes)
#plot_ordering_genes(cds)

## Trajectory Step 2: reducing the dimensionality of the data
cds <- reduceDimension(cds, max_components = 2,
                            method = 'DDRTree')

## Trajectory step 3: order cells along the trajectory
cds <- orderCells(cds)

GM_state <- function(cds){
  if (length(unique(pData(cds)$State)) > 1){
    T0_counts <- table(pData(cds)$State, pData(cds)$stage)[,"Early"]
    return(as.numeric(names(T0_counts)[which
                                       (T0_counts == max(T0_counts))]))
  } else {
    return (1)
  }
}
cds <- orderCells(cds, root_state = GM_state(cds))

## visualization of trajectory
plot_cell_trajectory(cds, color_by = "stage")

plot_cell_trajectory(cds, color_by = "State")
plot_cell_trajectory(cds, color_by = "State") +
  facet_wrap(~State, nrow = 1)

plot_cell_trajectory(cds, color_by = "Pseudotime")


## visualize gene expression
blast_genes <- row.names(subset(fData(cds),
                                gene_short_name %in% c("Il17a", 'Ifng')))
plot_genes_jitter(cds[blast_genes,],
                  grouping = "State",
                  min_expr = 0.1)


HSMM_expressed_genes <-  row.names(subset(fData(cds),
                                          num_cells_expressed >= 10))
HSMM_filtered <- cds[HSMM_expressed_genes,]
my_genes <- row.names(subset(fData(HSMM_filtered),
                             gene_short_name %in% c("Il17a", 'Ifng')))
cds_subset <- HSMM_filtered[my_genes,]
plot_genes_in_pseudotime(cds_subset, color_by = "stage")



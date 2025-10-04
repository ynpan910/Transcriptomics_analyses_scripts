#-------------------------------------------------monocle3----------------------------#
library(monocle3)
library(Seurat);
library(SeuratWrappers)
library(tidyverse)
source('color.R')


ab<- readRDS('seu.escape.rds')
metadata<- ab@meta.data
DefaultAssay(ab)


#Idents(ab)<-'cell_type_me'
#Idents(ab)<- factor(Idents(ab), levels = val1)
#Levels: HM HLA CRM IRM DAM Proliferating cell


#'----------------------------------------
#' from seurat to cds
#'----------------------------------------

cds<- as.cell_data_set(ab)
cds <- estimate_size_factors(cds) # this is super important! otherwise the y-axis outputed by plot_gene_pseudotime will be divided by 1*10-04
colData(cds)
fData(cds)
rownames(fData(cds))[1:10]
fData(cds)$gene_short_name<- rownames(fData(cds))


#'----------------------------------------
#' transfer clustring info from seurat's umap
#'----------------------------------------


#--------- transfer seurat_cluster info---------
# assign cluster info
cds@clusters$UMAP$clusters<- ab$seurat_clusters

# assign UMAP coordinate-cell embeddings
cds@int_colData@listData$reducedDims$UMAP<- ab@reductions$umap.cca@cell.embeddings


# ---------- plot cluster assignment ----------
cluster.before.trajectory<- plot_cells(cds,
                                       color_cells_by = 'cluster',
                                       label_groups_by_cluster = F,
                                       group_label_size = 5)+
  theme(legend.position = 'right')

cluster.names<- plot_cells(cds,
                           color_cells_by = 'cell_type_me',
                           label_groups_by_cluster = F,
                           group_label_size = 5)+
  theme(legend.position = 'right')

cluster.before.trajectory | cluster.names


# ---------- assign partitions----------

# give all cells as partition 1 first
recreate.partition<- c(rep(1, length(cds@colData@rownames)))
names(recreate.partition)<- cds@colData@rownames

# CRM on partition 2
#recreate.partition[names(recreate.partition) %in% rownames(metadata %>% filter(cell_type_me == 'CRM'))] <- 2


# IRM on partition 3
#recreate.partition[names(recreate.partition) %in% rownames(metadata %>% filter(cell_type_me == 'IRM'))] <- 3


# double check the partition assign
recreate.partition<- as.factor(recreate.partition)
table(recreate.partition)
nrow(metadata %>% filter(cell_type_me=='CRM'))
nrow(metadata %>% filter(cell_type_me=='IRM'))
nrow(metadata %>% filter(cell_type_me=='Proliferating cell'))


cds@clusters$UMAP$partitions<- recreate.partition
names(cds@clusters$UMAP$partitions)<- colnames(cds)





#'----------------------------------------
#' learn trajectory graph
#'----------------------------------------

cds<- learn_graph(cds, use_partition = T)

plot_cells(cds,
           color_cells_by = 'cell_type_me',
           group_cells_by= 'partition',
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F,
           group_label_size = 5)


# a helper function to identify the root principal points:
get_earliest_principal_node <- function(cds, cell_type_me="HM"){
  cell_ids <- which(colData(cds)[, "cell_type_me"] == cell_type_me)
  
  closest_vertex <-
    cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.matrix(closest_vertex[colnames(cds), ])
  root_pr_nodes <-
    igraph::V(principal_graph(cds)[["UMAP"]])$name[as.numeric(names
                                                              (which.max(table(closest_vertex[cell_ids,]))))]
  
  root_pr_nodes
}

get_earliest_principal_node(cds) #' "Y_51"
cds <- order_cells(cds, root_pr_nodes=get_earliest_principal_node(cds))
#' OR just simply choose the cluster but not using the programming method by : root_cells = colnames(cds[, colData(cds)$cell_type_me=='HM']))


# // trajectory analysis is now run done //


#'----------------------------------------
#' rough visualization
#'----------------------------------------

# ------------- plot UMAP colored by pseudotime -------------
plot_cells(cds,
           color_cells_by = 'cell_type_me',
           group_cells_by= 'cell_type_me',
           label_groups_by_cluster = T,
           label_branch_points = F,
           label_roots = T,
           label_leaves = F,
           group_label_size = 5,
           show_trajectory_graph = T,
           cell_size = 0.5)

plot_cells(cds,
           color_cells_by = 'pseudotime',
           group_cells_by= 'cell_type_me',
           label_groups_by_cluster = T,
           label_branch_points = F,
           label_roots = T,
           label_leaves = F,
           group_label_size = 5,
           show_trajectory_graph = T,
           cell_size = 0.5)

# ------------- cells ordered by monocle3 pseudotime using boxplot --------------#
pseudotime(cds)
cds$monocle3_pseudotime<- pseudotime(cds)
data.pseudo<- as.data.frame(colData(cds))

data.pseudo$monocle3_pseudotime_plot <- ifelse(
  is.infinite(data.pseudo$monocle3_pseudotime),
  max(data.pseudo$monocle3_pseudotime[is.finite(data.pseudo$monocle3_pseudotime)], na.rm = TRUE) + 1,
  data.pseudo$monocle3_pseudotime
)
pdf('cell_types_by_pseudotime_boxplot.pdf', width = 8, height = 5)
ggplot(data.pseudo,
       aes(x = reorder(cell_type_me, monocle3_pseudotime_plot, median),
           y = monocle3_pseudotime_plot,
           fill = cell_type_me)) +
  geom_boxplot() +
  xlab("") +
  ylab("Pseudotime")+
  theme_classic()+
  theme(axis.text = element_text(size = 14, face = 'bold', colour = 'black'),
        axis.title.y = element_text(size = 14, face = 'bold', colour = 'black'),
        legend.title = element_text(size = 14, face = 'bold', colour = 'black'),
        legend.text = element_text(size = 12, face = 'bold', colour = 'black'))+
  labs(fill = 'Cell type')
graphics.off()




# ------------- show ig gene expressions change as a function of pseudotime -------------#



markers<- c( 'TMEM119', 'CD9')

for (gene in markers) {
  p<- plot_genes_in_pseudotime(cds[gene,],color_cells_by="cell_type_me", min_expr = 0)+ylim(0, 3)
  #pdf(paste('ab', gene, 'overall_trajectory.pdf', sep = '_'))
  print(p)
  #graphics.off()
}

saveRDS(cds, 'ab_trajectory.rds')

              

       
                     
#'-----------------------------
#' ****visualization for pp ******                  
#'-----------------------------                    
library(monocle3)
library(Seurat);
library(SeuratWrappers)
library(tidyverse)
source('color.R')


# --------- a customized function for plotting trajectory of genes with pseudotime -------------#
source('monocle3_plot_gene_in_peudotime_zero.R')
#' this customized function is to plot a 'cleaner' version of the original plot_gene_in_pseudotime function. 
#' It keep the spline trajectory line but put all scattered points (cells) to be on y=0 axis. Therefore, the points (cells) on y=0 represents 'cell states' or 'cell tyeps'.

plot_genes_in_pseudotime_zero <- function (cds_subset, min_expr = NULL, cell_size = 0.75, nrow = NULL, 
                                           ncol = 1, panel_order = NULL, color_cells_by = "pseudotime", 
                                           trend_formula = "~ splines::ns(pseudotime, df=3)", label_by_short_name = TRUE, 
                                           vertical_jitter = NULL, horizontal_jitter = NULL) 
{

  colData(cds_subset)$pseudotime <- pseudotime(cds_subset)
  cds_subset = cds_subset[, is.finite(colData(cds_subset)$pseudotime)]
  cds_exprs <- SingleCellExperiment::counts(cds_subset)
  cds_exprs <- Matrix::t(Matrix::t(cds_exprs)/size_factors(cds_subset))
  cds_exprs <- reshape2::melt(round(as.matrix(cds_exprs)))
  
  colnames(cds_exprs) <- c("f_id", "Cell", "expression")
  cds_colData <- as.data.frame(colData(cds_subset))
  cds_rowData <- as.data.frame(rowData(cds_subset))
  cds_exprs <- merge(cds_exprs, cds_rowData, by.x = "f_id", by.y = "row.names")
  cds_exprs <- merge(cds_exprs, cds_colData, by.x = "Cell", by.y = "row.names")
  cds_exprs$adjusted_expression <- cds_exprs$expression
  cds_exprs$feature_label <- if (label_by_short_name && "gene_short_name" %in% names(cds_exprs)) {
    ifelse(is.na(cds_exprs$gene_short_name), cds_exprs$f_id, as.character(cds_exprs$gene_short_name))
  } else {
    cds_exprs$f_id
  }
  cds_exprs$feature_label <- factor(cds_exprs$feature_label)
  
  # fit trajectory
  new_data <- as.data.frame(colData(cds_subset))
  new_data$Size_Factor = 1
  model_tbl = fit_models(cds_subset, model_formula_str = trend_formula)
  model_expectation <- model_predictions(model_tbl, new_data = new_data)
  colnames(model_expectation) <- colnames(cds_subset)
  expectation <- plyr::ddply(cds_exprs, plyr::.(f_id, Cell), 
                             function(x) {
                               data.frame(expectation = model_expectation[x$f_id, x$Cell])
                             })
  cds_exprs <- merge(cds_exprs, expectation)
  cds_exprs$expectation[cds_exprs$expectation < min_expr] <- min_expr
  
  # --- plot ---
  q <- ggplot(aes(pseudotime), data = cds_exprs)
  if (!is.null(color_cells_by)) {
    # force all points to y=0
    q <- q + geom_point(aes_string(y = 0, color = color_cells_by),
                        size = I(cell_size),
                        position = position_jitter(horizontal_jitter, vertical_jitter))
    if (methods::is(colData(cds_subset)[, color_cells_by], "numeric")) {
      q <- q + viridis::scale_color_viridis(option = "C")
    }
  } else {
    q <- q + geom_point(aes(y = 0),
                        size = I(cell_size),
                        position = position_jitter(horizontal_jitter, vertical_jitter))
  }
  
  # keep spline trajectory line
  q <- q + geom_line(aes(x = pseudotime, y = expectation), data = cds_exprs)
  
  # axis/scales
  q <- q + facet_wrap(~feature_label, nrow = nrow, ncol = ncol, scales = "free_y")
  q <- q + ylab("Expression") + xlab("pseudotime")
  q <- q + theme_classic()
  q
}




# ---------- gene trajectory ----------
cds<- readRDS('ab_trajectory.rds')  
markers<- c( 'TMEM119', 'CD9', "HLA-DQB1", 'CCL4', 'IFITM3','MKI67')

for (gene in markers) {
  p<- plot_genes_in_pseudotime_zero(cds[gene,],color_cells_by="cell_type_me", min_expr = 0, horizontal_jitter = 0, vertical_jitter = 0)+
    scale_color_manual(values = col1)+
    labs(color = 'Cell type')+
    theme(legend.title  = element_text(face = 'bold', size = 14),
          legend.text = element_text(face = 'bold', size = 13),
          axis.title = element_text(face = 'bold', size = 13, colour = 'black'),
          axis.text = element_text(face = 'bold', size = 13, colour = 'black'),
          strip.text = element_text(face = 'bold', size = 14)  )+
    guides(color = guide_legend(override.aes = list(size = 5)))
    
  pdf(paste('ab', gene, 'overall_trajectory.pdf', sep = '_'), width = 6, height = 4)
  print(p)
  graphics.off()
}


# ---------- UMAP colored by pseudotime --------------
pdf('ab_pseudotime_UMAP_by_pseudotime_w_lines.pdf', width = 6, height = 4)
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = T,
           label_branch_points = F,
           label_roots = T,
           label_leaves = F,
           group_label_size = 5,
           show_trajectory_graph = T,
           cell_size = 1)
graphics.off()

pdf('ab_pseudotime_UMAP_by_pseudotime_wo_lines.pdf', width = 6, height = 4)
plot_cells(cds,
           color_cells_by = 'pseudotime',
           label_groups_by_cluster = T,
           label_branch_points = F,
           label_roots = F,
           label_leaves = F,
           group_label_size = 5,
           show_trajectory_graph = F,
           cell_size = 1)
graphics.off()


# ---------- UMAP colored by cell type --------------
pdf('ab_pseudotime_UMAP_by_cell_type_w_lines.pdf', width = 6, height = 4)
plot_cells(cds,
           color_cells_by = 'cell_type_me',
           #group_cells_by= 'cell_type_me',
           label_groups_by_cluster = F,
           label_branch_points = F,
           label_cell_groups = F,
           label_roots = F,
           label_leaves = F,
           group_label_size = 5,
           show_trajectory_graph = T,
           cell_size = 1)+
  scale_color_manual(values = col1)
graphics.off()


# ------------ horizontal barchart split by genotype ---------------
library(ggpubr) 
pseudotime(cds)
cds$monocle3_pseudotime<- pseudotime(cds)
data.pseudo<- as.data.frame(colData(cds))




# l1
l1<- data.pseudo %>% filter(cell_type_me %in% c('HM',  'HLA', 'CRM'))
l1$genotype<- factor(l1$genotype, levels = c('KO', 'WT'))
p<- ggplot(l1, aes(x = genotype, y = monocle3_pseudotime)) +
  geom_jitter(aes(color = cell_type_me), width = 0.2, size = 1.5, alpha =1) +  
  geom_boxplot(outlier.shape = NA, alpha = 0, width = 0.5, size = 1) + 
  labs(x = "Genotype", y = "Pseudotime", color = "Cell type") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 12),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 13, face = 'bold', color = 'black'),
    legend.position = "right"
  )+
  coord_flip()+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(values=c( 'HM' = '#1E79B4', 
                                'HLA' =  '#FF7F0B',
                               'CRM' =  '#DF156C',
                              'WT'= '#54B07C', 
                              'KO'= '#EE8432'))+
  stat_compare_means(
    comparisons = list(
      c("WT", "KO")),
    method = "t.test",
    label = "p.format"  # or "p.format" to show actual p-value
  )
pdf('ipsc_trajectory_box_by_genotype_l1.pdf', width = 6,height = 4)
p
graphics.off()






# l2
l2<- data.pseudo %>% filter(cell_type_me %in% c('HM',  'IRM'))
l2$genotype<- factor(l2$genotype, levels = c('KO', 'WT'))
p<- ggplot(l2, aes(x = genotype, y = monocle3_pseudotime)) +
  geom_jitter(aes(color = cell_type_me), width = 0.2, size = 1.5, alpha =1) +  
  geom_boxplot(outlier.shape = NA, alpha = 0, width = 0.5, size = 1) + 
  labs(x = "Genotype", y = "Pseudotime", color = "Cell type") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 12),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 13, face = 'bold', color = 'black'),
    legend.position = "right"
  )+
  coord_flip()+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(values=c( 'HM' = '#1E79B4', 
                               'IRM' =  '#D72626',
                               'WT'= '#54B07C', 
                               'KO'= '#EE8432'))+
  stat_compare_means(
    comparisons = list(
      c("WT", "KO")),
    method = "t.test",
    label = "p.format"  # or "p.format" to show actual p-value
  )
pdf('ipsc_trajectory_box_by_genotype_l2.pdf', width = 6,height = 4)
p
graphics.off()




# l3
l3<- data.pseudo %>% filter(cell_type_me %in% c('HM',  'DAM'))
l3$genotype<- factor(l3$genotype, levels = c('KO', 'WT'))
p<- ggplot(l3, aes(x = genotype, y = monocle3_pseudotime)) +
  geom_jitter(aes(color = cell_type_me), width = 0.2, size = 1.5, alpha =1) +  
  geom_boxplot(outlier.shape = NA, alpha = 0, width = 0.5, size = 1) + 
  labs(x = "Genotype", y = "Pseudotime", color = "Cell type") +
  theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size = 12),
    axis.title.x = element_text(size = 14, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 14, face = "bold", color = "black"),
    legend.title = element_text(size = 14, face = "bold", color = "black"),
    legend.text = element_text(size = 13, face = 'bold', color = 'black'),
    legend.position = "right"
  )+
  coord_flip()+
  guides(color = guide_legend(override.aes = list(size = 5)))+
  scale_color_manual(values=c( 'HM' = '#1E79B4', 
                               'DAM' =   '#149759',
                               'WT'= '#54B07C', 
                               'KO'= '#EE8432'))+
  stat_compare_means(
    comparisons = list(
      c("WT", "KO")),
    method = "t.test",
    label = "p.format"  # or "p.format" to show actual p-value
  )
pdf('ipsc_trajectory_box_by_genotype_l3.pdf', width = 6,height = 4)
p
graphics.off()



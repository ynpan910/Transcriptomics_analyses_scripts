
#---------------------------------make seurat objects for each sample--------------------------------#

library(patchwork)
library(DropletUtils)
library(Seurat)


#---repeat the following for each sample----#

filter_matrix <- Read10X("./Sample_3792-FM-S2-A/filtered_feature_bc_matrix/")

write10xCounts("./Sample_3792-FM-S2-A/filtered_feature_bc_matrix.h5", filter_matrix, type = "HDF5",
               genome = "mm10", version = "3", overwrite = TRUE,
               gene.id = rownames(filter_matrix),
               gene.symbol = rownames(filter_matrix))

#load in
m11<- Load10X_Spatial(data.dir = './Sample_3792-FM-S2-A/', filename = 'filtered_feature_bc_matrix.h5')

metadata_m11<- m11@meta.data

plot1 <- VlnPlot(m11, features = "nCount_Spatial", pt.size = 0.1) + NoLegend()
plot2 <- SpatialFeaturePlot(m11, features = "nCount_Spatial") + theme(legend.position = "right")
wrap_plots(plot1, plot2)


m11 <- SCTransform(m11, assay = "Spatial", verbose = FALSE)

m11 <- RunPCA(m11, assay = "SCT", verbose = FALSE)
m11 <- FindNeighbors(m11, reduction = "pca", dims = 1:30)
m11 <- FindClusters(m11, verbose = FALSE)
m11 <- RunUMAP(m11, reduction = "pca", dims = 1:30)

p1 <- DimPlot(m11, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(m11, label = TRUE, label.size = 3)
p1 + p2

metadata_m11<- m11@meta.data

saveRDS(m11,file='Male_rep1_sec1_spatial_seurat_object.RDS')



#--------------------------------Deconvolution to see cell composition of each spot-----------------#


library(tidyverse)
library(Seurat)
library(Giotto)

my_python_path<- "C:/Users/yining.pan/miniconda3/python.exe"
instr<- createGiottoInstructions(python_path = my_python_path)


#---repeat the following for each sample----#


#--------------------
# process st data |
#--------------------
#' convert seurat object to Giotto object

m11<- readRDS('Male_rep1_sec1_spatial_seurat_object.RDS')

spatial.exp <- m11[['Spatial']]$counts
spatial.loc <- GetTissueCoordinates(m11)
rownames(spatial.loc)<-NULL
spatial.loc<- spatial.loc %>% column_to_rownames(var = 'cell')
identical(colnames(spatial.exp),rownames(spatial.loc))

## Generate Giotto objects and cluster spots
library(Giotto)
m11.giotto <- createGiottoObject(
  expression  = spatial.exp,
  spatial_locs = spatial.loc, instructions = instr)

m11.giotto <- filterGiotto(
  gobject = m11.giotto,
  expression_threshold = 0,
  feat_det_in_min_cells = 0,
  min_det_feats_per_cell = 0,
  expression_values = c('raw'),
  verbose = T)

m11.giotto <- normalizeGiotto(gobject = m11.giotto)
m11.giotto <- calculateHVF(gobject = m11.giotto)

gene_metadata = fDataDT(m11.giotto)

#featgenes = gene_metadata[hvf == 'yes']$gene_ID

m11.giotto <- runPCA(
  gobject = m11.giotto,
  scale_unit = F)

signPCA(m11.giotto,  scale_unit = F)

m11.giotto <- createNearestNetwork(gobject = m11.giotto, dimensions_to_use = 1:10, k = 10)

m11.giotto <- doLeidenCluster(gobject = m11.giotto, resolution = 0.4, n_iterations = 1000,python_path = "C:/Users/yining.pan/miniconda3/python.exe")

saveRDS(m11.giotto, 'm11.giotto.rds')


#-------------------
# spatialDWLS run |
#-------------------
#' takes five hours to finish

m11.giotto<- readRDS('m11.giotto.rds')
Sig_exp<- readRDS('Sig_exp_spatialDWLS_malemouse.RDS')


library(Rfast)
library(quadprog)
m11.decon <- runDWLSDeconv(
  gobject = m11.giotto, sign_matrix = Sig_exp)


saveRDS(m11.decon, 'spatialDWLS_deconvoluted_m11.RDS')


#----------------------------------------------
# Visualization of each deconvoluted sample |
#----------------------------------------------

## bar chart
m11.decon<- readRDS('spatialDWLS_deconvoluted_m11.RDS')
decon.res.m11<- as.data.frame(m11.decon@spatial_enrichment$cell$rna$DWLS@enrichDT)

fo1 <- colMeans(decon.res.m11[, 2:(ncol(decon.res.m11) - 1)])
fo2 <- apply(decon.res.m11[, 2:(ncol(decon.res.m11) - 1)], 2, sd)
fo3 <- rbind(fo1, fo2)
rownames(fo3) <- c('mean', 'sd')
fo4 <- data.frame(t(fo3))
fo4$cell <- rownames(fo4)
rownames(fo4) <- NULL
fo4


fo4<- fo4 %>% arrange(desc(mean))
fo4$cell<- factor(fo4$cell, levels = unique(fo4$cell))

ggplot(fo4) +
  geom_bar(aes(x=cell, y=mean), stat="identity",
           width=0.5,
           position=position_dodge()) +
  theme_light() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))+
  labs(x='cell population', y='mean (%)')

## pie chart

library(Giotto)
library(ggplot2)
library(scatterpie)
library(data.table)


m11.giotto<- readRDS('m11.giotto.rds')
m11.decon<- readRDS('spatialDWLS_deconvoluted_m11.RDS')


plot_data<- as.data.frame(m11.decon@spatial_enrichment$cell$rna$DWLS@enrichDT)
library(tidyverse)
plot_data<- plot_data %>% column_to_rownames(var='cell_ID')


plot_data$`Excitatory neuron` <- rowSums(plot_data[, c(1:6, 9, 17, 19, 25, 26,  31:35, 37:39)], na.rm = TRUE)
plot_data$`Inhibitory neuron`<- rowSums( plot_data[, c(7, 8, 10, 11, 13, 14, 23)], na.rm = TRUE)
plot_data$Others<- rowSums(plot_data[, c(15, 22, 23, 27:30, 36, 40:42)], na.rm = TRUE)

names(plot_data)[names(plot_data) == 'Oligo'] <- 'Oligodendrocyte'
names(plot_data)[names(plot_data) == 'Astro'] <- 'Astrocyte'
names(plot_data)[names(plot_data) == 'Endo'] <- 'Endothelial cell'

plot_col <- c('Excitatory neuron', 'Inhibitory neuron', 'Oligodendrocyte', 'Astrocyte', 'Micro-PVM', 'Endothelial cell',  
              'VLMC', 'SMC-Peri', 'Others')

plot_data$x <- as.numeric(as.character(m11.giotto@spatial_locs$cell$raw$sdimx))
plot_data$y <- as.numeric(as.character(m11.giotto@spatial_locs$cell$raw$sdimy))



df <- data.frame()
p <- ggplot(df) + geom_point() + xlim(min(plot_data$x)-1, max(plot_data$x)+1) + ylim(min(plot_data$y)-1, max(plot_data$y)+1)


library(ggsci)
pdf('m11.visium.spatialDWLS.decon.pie.pdf', width = 6, height = 4)

p + geom_scatterpie(aes(x=x, y=y), data=plot_data, cols=plot_col, alpha=1, pie_scale = 0.4, color=NA) +theme_classic()+
  scale_fill_observable(name='Cell type')+theme(
    axis.line = element_blank(),        # Remove axis lines
    axis.text.x = element_blank(),         # Remove x-axis text
    axis.text.y = element_blank(),         # Remove y-axis text
    axis.title.x = element_blank(),        # Remove x-axis title
    axis.title.y = element_blank(),        # Remove y-axis title
    axis.ticks = element_blank(), 
    legend.title = element_text(size = 12, face = "bold", color = "black"),
    legend.text = element_text(size = 12, face = 'bold', color = 'black'))

graphics.off()


#----------------------------
#comparing cell composition |
#----------------------------

df<- read.csv('s222_st_public datasets.csv')

df$sex<- factor(df$sex, levels = c('Male', 'Female'))

library(tidyverse)
ggplot(df, aes(x=sex, y=Freq,fill=sex))+geom_boxplot(outlier.shape = NA) + geom_jitter(width = 0.2, color = "black")+
  scale_fill_manual(values = c('#0072B5FF', '#E18727FF')) +  theme_classic()+
  theme(
    axis.text.x = element_text(face = "bold", color = "black", size=13 ),
    axis.title.x = element_text(size = 12, face = "bold", color = "black"),
    axis.text.y = element_text(size = 12, face = "bold", color = "black"),
    axis.title.y = element_text(size = 12, face = "bold", color = "black"),
    legend.title = element_text(size = 13, face = "bold", color = "black"),
    legend.text = element_text(size = 13, face = 'bold', color = 'black'),
    legend.position = 'none',
    strip.text = element_text(face = "bold", size=12))+labs(x='', y='Oligodendrocyte composition (%)')


res.ftest <- var.test(Freq~sex, data = df)
res.ftest
t.res <- t.test(Freq ~ sex, data = df, var.equal = T)
t.res


#----------------------------------------------DEG and Enrichment analysis-------------------------#

m11<- readRDS('Male_rep1_sec1_spatial_seurat_object.RDS')
m12<- readRDS('Male_rep1_sec2_spatial_seurat_object.RDS')
m21<- readRDS('Male_rep2_sec1_spatial_seurat_object.RDS')
m22<- readRDS('Male_rep2_sec2_spatial_seurat_object.RDS')

f3<- readRDS('Female_3_spatial_seurat_object.RDS')
f4<- readRDS('Female_4_spatial_seurat_object.RDS')
f5<- readRDS('Female_5_spatial_seurat_object.RDS')

ls<- list(m11, m12 ,m21, m22, f3, f4, f5)

ls.name<- c('m11', 'm12' ,'m21', 'm22', 'f3', 'f4', 'f5')

brain.merge<- merge(ls[[1]], y=ls[2:7],
                    add.cell.ids = ls.name)

DefaultAssay(brain.merge) <- "SCT"
VariableFeatures(brain.merge) <- c(VariableFeatures(m11), VariableFeatures(m12), 
                                   VariableFeatures(m21), VariableFeatures(m22),
                                   VariableFeatures(f3), VariableFeatures(f4),
                                   VariableFeatures(f5))
brain.merge <- RunPCA(brain.merge, verbose = FALSE)
brain.merge <- FindNeighbors(brain.merge, dims = 1:30)
brain.merge <- FindClusters(brain.merge, verbose = FALSE, resolution = 0.5)
brain.merge <- RunUMAP(brain.merge, dims = 1:30)

DimPlot(brain.merge, reduction = "umap", label = TRUE)
SpatialDimPlot(brain.merge, label = T, label.size = 3, images = 'slice1')
SpatialDimPlot(brain.merge, label = T, label.size = 3, images = 'slice1', image.alpha = 0)


#--------------------
# spot annotation
#--------------------
my.dat<- read.csv('~/AD_proposal/cell_type_marker.csv', stringsAsFactors = F, row.names = NULL, check.names = F)
my.mark<- apply(my.dat, 2, function(x) x[x != ''])



library(stringr)
library(tidyverse)
pdf('feature_heatmap1.pdf', width = 10, height = 4)
for (i in 1:length(my.mark)) {
  heatmap_plot <- DoHeatmap(brain.merge, features = str_to_title(my.mark[[i]])) +
    ggtitle(names(my.mark)[i]) 
  print(heatmap_plot)
}
graphics.off()


#-------------------------
# cell_type_a area plot
#--------------------------
library(RColorBrewer)
library(tidyverse)
mycolors = c(brewer.pal(name="Dark2", n = 8), brewer.pal(name="Paired", n = 8))

## male
p1<- SpatialDimPlot(brain.merge, label = T, label.size = 3, images = 'slice1.2', image.alpha = 0, pt.size.factor = 4) & scale_fill_manual(values = mycolors)


p2<- SpatialFeaturePlot(brain.merge, features = c("Gad1", "Gad2", 'Lhx6'),images = 'slice1.2',image.alpha = 0, pt.size.factor = 4)

pdf('Male_1_2_cell_type_a_expression.pdf', width = 12, height = 5)
p1+p2
graphics.off()

## female

p1<- SpatialDimPlot(brain.merge, label = T, label.size = 3, images = 'slice1.6', image.alpha = 0, pt.size.factor = 4) & scale_fill_manual(values = mycolors)


p2<- SpatialFeaturePlot(brain.merge, features = c("Gad1", "Gad2", 'Lhx6'),images = 'slice1.6',image.alpha = 0, pt.size.factor = 4)

pdf('Female_4_cell_type_a_expression.pdf', width = 12, height = 5)
p1+p2
graphics.off()

#----------------------------------
# cell_type_a DEG AD Female vs. AD male
#-----------------------------------
brain.merge<- RenameIdents(brain.merge, 
                           '1'='Cell_type_a', 
                           '13'='Cell_type_a')
brain.merge$cell_type_me<- brain.merge@active.ident
brain.merge.cell_type_a<- subset(brain.merge, subset = cell_type_me=='Cell_type_a')

metadata.cell_type_a<- brain.merge.cell_type_a@meta.data
metadata.cell_type_a$sampleID<- sub("_.*", "", rownames(metadata.cell_type_a))
metadata.cell_type_a$sex<- substr(metadata.cell_type_a$sampleID, 1, 1)
metadata.cell_type_a$sex<- ifelse(metadata.cell_type_a$sex=='m', 'Male', 'Female')
metadata.cell_type_a$sex<- factor(metadata.cell_type_a$sex, levels = c('Male', 'Female'))
brain.merge.cell_type_a@meta.data<- metadata.cell_type_a

## pseudobulkDEG

brain.merge.cell_type_a[['Spatial']]<- JoinLayers(brain.merge.cell_type_a[['Spatial']])

### saveRDS(brain.merge.cell_type_a, 'brain.merge.cell_type_a.rds')


rawcount<- brain.merge.cell_type_a[['Spatial']]$counts
metadata<- brain.merge.cell_type_a@meta.data

library(scran)
sce<- SingleCellExperiment(assays=list(counts=rawcount), colData=metadata)
sce
groups<- colData(sce)[, c('cell_type_me', 'sampleID')]
summed<- aggregateAcrossCells(sce, groups)
summed

#colData(summed)$genotype<- as.factor(colData(summed)$genotype)
#colData(summed)$sex<- as.factor(colData(summed)$sex)
#colData(summed)$cell_type_me<- as.factor(colData(summed)$cell_type_me)

de.results<- pseudoBulkDGE(summed, label=summed$cell_type_me, design=~sex+ncells, condition=summed$sex, coef='sexFemale')

deg.cella<- as.data.frame(de.results[['Cell_type_a']])

deg.cella<- deg.cella[order(deg.cella$PValue),]
deg.cella$ratio<- 2^deg.cella$logFC
deg.cella$fold_change<- ifelse(deg.cella$ratio>1, deg.cella$ratio, -1/deg.cella$ratio)

###write.csv(deg.cella, 'deg.cella.FvsM.AD.csv')

deg.cella.sig<- deg.cella %>% filter(FDR<0.05 & abs(fold_change)>=1.2)
deg.cella.sig.up<- deg.cella.sig %>% filter(fold_change>0)
deg.cella.sig.dn<- deg.cella.sig %>% filter(fold_change<0)


#-------------------------
# ORA on DEG
#-------------------------
library(clusterProfiler)
library(org.Mm.eg.db)

dn.gene<- bitr(rownames(deg.cella.sig.dn), fromType = "SYMBOL",
               toType = "ENTREZID",
               OrgDb = org.Mm.eg.db)
up.gene<- bitr(rownames(deg.cella.sig.up), fromType = "SYMBOL", toType = "ENTREZID",OrgDb = org.Mm.eg.db)

## GO
up.go <- enrichGO(
  gene           = up.gene$ENTREZID,
  OrgDb          = org.Mm.eg.db,
  ont            = "BP",
  pAdjustMethod  = "BH"
)
dn.go <- enrichGO(
  gene           = dn.gene$ENTREZID,
  OrgDb          = org.Mm.eg.db,
  ont            = "BP",
  pAdjustMethod  = "BH"
)

up.go.res<- up.go@result
dn.go.res<- dn.go@result

write.csv(up.go.res, file=paste0('up.go.', 'Cella', '.csv'))
write.csv(dn.go.res, file=paste0('dn.go.', 'Cella', '.csv'))

## Reactome
library(ReactomePA)
library(tidyverse)
library(clusterProfiler)
library(org.Mm.eg.db)

dn.reactome<- enrichPathway(gene=dn.gene$ENTREZID, organism = 'mouse')
up.reactome<- enrichPathway(gene=up.gene$ENTREZID, organism = 'mouse')

dn.reactome.res<- dn.reactome@result
up.reactome.res<- up.reactome@result

write.csv(up.reactome.res, file=paste0('up.reactome.', 'Cella', '.csv'))
write.csv(dn.reactome.res, file=paste0('dn.reactome.', 'Cella', '.csv'))

## KEGG
up.kegg <- enrichKEGG(gene = up.gene$ENTREZID, organism = 'mmu')
dn.kegg <- enrichKEGG(gene = dn.gene$ENTREZID, organism = 'mmu')

dn.kegg.res<- dn.kegg@result
up.kegg.res<- up.kegg@result

write.csv(up.kegg.res, file=paste0('up.kegg.', 'Cella', '.csv'))
write.csv(dn.kegg.res, file=paste0('dn.kegg.', 'Cella', '.csv'))


#------------------
# visualization
#------------------

library(patchwork)
library(stringr)
library(tidyverse)

## GO

up.go.res<- read.csv('visium_AD_FvsM_up.go.Cella.csv')
dn.go.res<- read.csv('visium_AD_FvsM_dn.go.Cella.csv')


top_positive <- up.go.res[order(up.go.res$p.adjust), ][1:15, ] 
top_positive$group<- 'up'
top_negative <- dn.go.res[order(dn.go.res$p.adjust), ][1:15, ]
top_negative$group<- 'dn'


df<- rbind(top_positive, top_negative)
df$Description<- str_to_sentence(df$Description)

df$`-log10padjust`<- -log10(df$p.adjust)
df$`-log10padjust_2`<- ifelse(df$group=='up',df$`-log10padjust`, -df$`-log10padjust`)

df$Description<- factor(df$Description, levels = rev(df$Description))

p<- ggplot(df,
           aes(x=`-log10padjust_2`, y=Description, fill=group))+
  geom_col()+
  theme_bw()
p


mytheme<- theme(
  legend.position = 'none',
  axis.text.y=element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(colour = 'grey60', size = 1.1),
  axis.text =  element_text(size = 12, face = 'bold', color = 'black')
  
)

p1<- p+mytheme
p1

upp<- df[which(df$`-log10padjust_2`>0),]
downn<- df[which(df$`-log10padjust_2`<0),]

p2<- p1+geom_text(
  data = upp,
  aes(x=-0.2, y=Description, label=Description),
  size=4,
  hjust=1,
  fontface = "bold" 
)
p2

p3<- p2+geom_text(
  data = downn,
  aes(x=0.2, y=Description, label=Description),
  size=4,
  hjust=0,
  fontface = "bold" 
)
p3


p4<- p3+labs(x='-Log10(padj)', y='', title = '')+
  theme(plot.title = element_text(hjust = 0.5, size = 14))+
  scale_fill_manual(values = c('#1C1BF5', '#9D1F14' ))+
  xlim(c(-7.5, 3))
p4

pdf('visium_AD_FvsM_cella_ORA_GO.pdf', width = 15, height = 8)
p4
graphics.off()



## Reactome

up.go.res<- read.csv('visium_AD_FvsM_up.reactome.Cella.csv')
dn.go.res<- read.csv('visium_AD_FvsM_dn.reactome.Cella.csv')


top_positive <- up.go.res[order(up.go.res$pvalue), ][1:15, ] 
top_positive$group<- 'up'
top_negative <- dn.go.res[order(dn.go.res$pvalue), ][1:15, ]
top_negative$group<- 'dn'


df<- rbind(top_positive, top_negative)
df$Description<- str_to_sentence(df$Description)

df$`-log10pvalue`<- -log10(df$pvalue)
df$`-log10pvalue_2`<- ifelse(df$group=='up',df$`-log10pvalue`, -df$`-log10pvalue`)

df$Description<- factor(df$Description, levels = rev(df$Description))

p<- ggplot(df,
           aes(x=`-log10pvalue_2`, y=Description, fill=group))+
  geom_col()+
  theme_bw()
p

mytheme<- theme(
  legend.position = 'none',
  axis.text.y=element_blank(),
  axis.ticks.y = element_blank(),
  panel.grid.major = element_blank(),
  panel.grid.minor = element_blank(),
  panel.border = element_blank(),
  axis.line.x = element_line(colour = 'grey60', size = 1.1),
  axis.text =  element_text(size = 12, face = 'bold', color = 'black')
  
)

p1<- p+mytheme
p1

upp<- df[which(df$`-log10pvalue_2`>0),]
downn<- df[which(df$`-log10pvalue_2`<0),]


p2<- p1+geom_text(
  data = upp,
  aes(x=-0.2, y=Description, label=Description),
  size=4,
  hjust=1,
  fontface = "bold" 
)
p2

p3<- p2+geom_text(
  data = downn,
  aes(x=0.2, y=Description, label=Description),
  size=4,
  hjust=0,
  fontface = "bold" 
)
p3


p4<- p3+labs(x='-Log10(pvalue)', y='', title = '')+
  theme(plot.title = element_text(hjust = 0.5, size = 14))+
  scale_fill_manual(values = c('#1C1BF5', '#9D1F14' ))+
  xlim(c(-4, 4))
p4

pdf('visium_AD_FvsM_cella_ORA_reactome.pdf', width = 15, height = 8)
p4
graphics.off()

















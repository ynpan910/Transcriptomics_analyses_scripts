#-------------------------------1. preprocessing steps in R---------------------------------#
library(SCENIC)
#packageVersion('SCENIC')
##have downloaded features datasets



setwd("~/AD_proposal/scenic_preprocessing")
all<- readRDS('~/AD_proposal/all.rds')
library(Seurat)
f<- subset(all, subset=msex=='Female')
m<- subset(all, subset=msex=='Male')

f_abcd_ad<- subset(f, subset=abcd_genotype %in% c('33', '34','44') & cogdx == '4')
m_abcd_ad<- subset(m, subset=abcd_genotype %in% c('33', '34','44') & cogdx == '4')

f_abcd_ad.metadata<- f_abcd_ad@meta.data
m_abcd_ad.metadata<- m_abcd_ad@meta.data

saveRDS(m_abcd_ad,'m_abcd_ad.rds')


#----repeat the following for male and female separately--------------#
cellInfo<- m_abcd_ad@meta.data
exprMat<- m_abcd_ad@assays$RNA$counts
dim(exprMat)
head(cellInfo)

identical(colnames(exprMat), rownames(cellInfo))

cellInfo <- data.frame(cellInfo)
cbind(table(cellInfo$cell_type_me))

#---------------------------------
# prepare the databases|
#---------------------------------

org <- "hgnc" # for human
dbDir <- "/home/yining.pan/AD_proposal" # RcisTarget databases location
myDatasetTitle <- "SCENIC example on Human brain" # choose a name for your analysis
data(defaultDbNames)
data(list="motifAnnotations_hgnc_v9", package="RcisTarget")

#rename the motif annnotion by attributing it to the variable that is in the error
motifAnnotations_hgnc <- motifAnnotations_hgnc_v9
dbs <- defaultDbNames[[org]]
scenicOptions <- initializeScenic(org=org, dbDir=dbDir, dbs=dbs, 
                                  datasetTitle=myDatasetTitle, nCores=10) 

saveRDS(cellInfo, file="int/cellInfo.Rds")
saveRDS(colVars, file="colVars.Rds")


scenicOptions@inputDatasetInfo$cellInfo <- "cellInfo.Rds"
scenicOptions@inputDatasetInfo$colVars <- "colVars.Rds"

#---------------------------------
# filtering low expressed genes |
#---------------------------------

## redefine the geneFiltering function
geneFiltering <- function(exprMat, scenicOptions,
                          minCountsPerGene=3*.01*ncol(exprMat),
                          minSamples=ncol(exprMat)*.01)
{
  # Load options: outFile_genesKept and dbFilePath
  outFile_genesKept <- NULL
  dbFilePath <- NULL
  if(class(scenicOptions) == "ScenicOptions")
  {
    dbFilePath <- getDatabases(scenicOptions)[[1]]
    outFile_genesKept <- getIntName(scenicOptions, "genesKept")
  }else{
    dbFilePath <- scenicOptions[["dbFilePath"]]
    outFile_genesKept <- scenicOptions[["outFile_genesKept"]]
  }
  if(is.null(dbFilePath)) stop("dbFilePath")
  
  # Check expression matrix (e.g. not factor)
  if(is.data.frame(exprMat)) 
  {
    supportedClasses <- paste(gsub("AUCell_buildRankings,", "", methods("AUCell_buildRankings")), collapse=", ")
    supportedClasses <- gsub("-method", "", supportedClasses)
    
    stop("'exprMat' should be one of the following classes: ", supportedClasses, 
         "(data.frames are not supported. Please, convert the expression matrix to one of these classes.)")
  }
  if(any(table(rownames(exprMat))>1))
    stop("The rownames (gene id/name) in the expression matrix should be unique.")
  
  # Calculate stats
  nCountsPerGene <- Matrix::rowSums(exprMat, na.rm = T)
  nCellsPerGene <- Matrix::rowSums(exprMat>0, na.rm = T)
  
  ## Show info
  message("Maximum value in the expression matrix: ", max(exprMat, na.rm=T))
  message("Ratio of detected vs non-detected: ", signif(sum(exprMat>0, na.rm=T) / sum(exprMat==0, na.rm=T), 2))
  message("Number of counts (in the dataset units) per gene:")
  print(summary(nCountsPerGene))
  message("Number of cells in which each gene is detected:")
  print(summary(nCellsPerGene))
  
  ## Filter
  message("\nNumber of genes left after applying the following filters (sequential):")
  # First filter
  # minCountsPerGene <- 3*.01*ncol(exprMat)
  genesLeft_minReads <- names(nCountsPerGene)[which(nCountsPerGene > minCountsPerGene)]
  message("\t", length(genesLeft_minReads), "\tgenes with counts per gene > ", minCountsPerGene)
  
  # Second filter
  # minSamples <- ncol(exprMat)*.01
  nCellsPerGene2 <- nCellsPerGene[genesLeft_minReads]
  genesLeft_minCells <- names(nCellsPerGene2)[which(nCellsPerGene2 > minSamples)]
  message("\t", length(genesLeft_minCells), "\tgenes detected in more than ",minSamples," cells")
  
  # Exclude genes missing from database:
  library(RcisTarget)
  motifRankings <- importRankings(dbFilePath) # either one, they should have the same genes
  genesInDatabase <- colnames(getRanking(motifRankings))
  
  genesLeft_minCells_inDatabases <- genesLeft_minCells[which(genesLeft_minCells %in% genesInDatabase)]
  message("\t", length(genesLeft_minCells_inDatabases), "\tgenes available in RcisTarget database")
  
  genesKept <- genesLeft_minCells_inDatabases
  if(!is.null(outFile_genesKept)){ 
    saveRDS(genesKept, file=outFile_genesKept)
    if(getSettings(scenicOptions, "verbose")) message("Gene list saved in ", outFile_genesKept)
  }
  return(genesKept)
}


genesKept <- geneFiltering(exprMat, scenicOptions=scenicOptions,
                           minCountsPerGene=25*.01*ncol(exprMat),
                           minSamples=ncol(exprMat)*.1)


interestingGenes <- c("ABCD", "TREM2")
# any missing?
interestingGenes[which(!interestingGenes %in% genesKept)]

exprMat_filtered <- exprMat[genesKept, ]
dim(exprMat_filtered)  #5382*28690

#------------------------------------------------------
# output to be imported in the next pyscenic pipeline |
#------------------------------------------------------

exprMat_filtered <- log2(exprMat_filtered+1) 
write.csv(exprMat_filtered, 'abcd_subset_m_exprMat_filtered_log2.csv')


#-------------------------------2. Generating sample.loom file in Python---------------------------------#
import loompy as lp
import pandas as pd
import scanpy as sc
import numpy as np

# load in csc file and check its basic information
x=pd.read_csv('abcd_subset_m_exprMat_filtered_log2.csv', header=0, index_col=0)
x.head()
x.shape

row_attrs={"Gene":np.array(x.index)}
row_attrs

col_attrs={"CellID": np.array(x.columns)}
col_attrs

x.values

# generate loom file
lp.create("sample.loom", x.values, row_attrs, col_attrs)

# read in and double check the loom file
import loompy

with loompy.connect("sample.loom") as ds:
  print("Row attributes:", ds.row_attrs.keys())  # Should include 'Gene'
  print("Column attributes:", ds.col_attrs.keys())  # Should include 'cellID'

  # Access attributes
  genes = ds.row_attrs["Gene"]
  cell_ids = ds.col_attrs["CellID"]
  print("Genes:", genes)
  print("Cell IDs:", cell_ids)


  
  
  
  
#-------------------------------3. TF inferring in Terminal---------------------------------#
## prerequisite: downloaded TF.txt, motif, and genes_vs_motifs
  

# step 1: GRN creating
  
cd '/home/yining.pan/scenic'
ls

arboreto_with_multiprocessing.py  --num_workers 20 --output adj.sample.tsv --method grnboost2 sample.loom hs_hgnc_tfs.txt
  

# step 2: cisTarget

tfs=./hs_hgnc_tfs.txt
feather=./hg19-500bp-upstream-7species.mc9nr.genes_vs_motifs.rankings.feather
tbl=./motifs-v9-nr.hgnc-m0.001-o0.0.tbl
input_loom=./sample.loom
ls

pyscenic ctx adj.sample.tsv $feather --annotations_fname $tbl --expression_mtx_fname $input_loom --mode “dask_multiprocessing” --output reg.csv --num_workers 20 --mask_dropouts

ls


# step 3: AUcell 

input_loom=./sample.loom
  
pyscenic aucell $input_loom reg.csv --output out_SCENIC.loom --num_workers 3
  
  
  
  
  
  
  
  
  

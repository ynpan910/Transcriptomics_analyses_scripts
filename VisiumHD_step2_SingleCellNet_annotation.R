library(singleCellNet)
# load query data
object<- readRDS('object.final.rds')

stQuery<- object@meta.data
expQuery<- object@assays$Spatial.008um$counts

genesQuery = rownames(expQuery)

rm(expQuery); rm(object)
gc()

#load reference data
load('Seurat.ss.rda')
expAllen <- ss.seurat@assays$RNA$counts

stAllen <- ss.seurat@meta.data

rm(ss.seurat)
stAllen<-droplevels(stAllen)


#Find genes in common to the data sets and limit analysis to these
commonGenes = intersect(rownames(expAllen), genesQuery)
length(commonGenes)
#[1] 17466
stAllen<- stAllen[!is.na(stAllen$subclass_label),]

#Split for training and assessment, and transform training data
set.seed(100) #can be any random seed number
stList = splitCommon(sampTab=stAllen, ncells=100, dLevel="subclass_label")
stTrain = stList[[1]]
expTrain = expAllen[commonGenes,rownames(stTrain)]

#train the classifier
stTrain$cell<- rownames(stTrain)
class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, nTopGenes = 10, nRand = 70, nTrees = 1000, nTopGenePairs = 25, dLevel = "subclass_label", colName_samp = "cell")

#Apply to held out data
##validate data
stTestList = splitCommon(sampTab=stList[[2]], ncells=100, dLevel="subclass_label") #normalize validation data so that the assessment is as fair as possible
stTest = stTestList[[1]]
expTest = expAllen[commonGenes,rownames(stTest)]

##predict
classRes_val_all = scn_predict(cnProc=class_info[['cnProc']], expDat=expTest, nrand = 50)

stTest$cell<- rownames(stTest)
tm_heldoutassessment = assess_comm(ct_scores = classRes_val_all, stTrain = stTrain, stQuery = stTest, dLevelSID = "cell", classTrain = "subclass_label", 
                                   classQuery = "subclass_label", 
                                   nRand = 50)

plot_PRs(tm_heldoutassessment)
plot_metrics(tm_heldoutassessment)


#Create a name vector label used later in classification heatmap where the values are cell types/ clusters and names are the sample names

nrand = 50
sla = as.vector(stTest$subclass_label)
names(sla) = as.vector(stTest$cell)
slaRand = rep("rand", nrand) 
names(slaRand) = paste("rand_", 1:nrand, sep='')
sla = append(sla, slaRand) #include in the random cells profile created

sc_hmClass(classMat = classRes_val_all,grps = sla, max=300, isBig=TRUE)

plot_attr(classRes=classRes_val_all, sampTab=stTest, nrand=nrand, dLevel="subclass_label", sid="cell")

# apply to query data

#objects_to_remove <- setdiff(ls(), "class_info")
#rm(list = objects_to_remove); gc()

rm(expAllen); gc()

## load query data
#object<- readRDS('object.final.rds')
#expQuery<- object@assays$Spatial.008um$counts
# to save memory, i saved this before hand and import here
expQuery<- readRDS('expQuery.rds')

expQuery<- expQuery[commonGenes,]

nqRand = 50
crParkall<-scn_predict(class_info[['cnProc']], expQuery, nrand=nqRand)
saveRDS(crParkall, 'crQueryall.rds')



stQuery$sample_name<- rownames(stQuery)
#visualization
sgrp = as.vector(stQuery$seurat_clusters)
names(sgrp) = as.vector(rownames(stQuery))
grpRand =rep("rand", nqRand)
names(grpRand) = paste("rand_", 1:nqRand, sep='')
sgrp = append(sgrp, grpRand)

# heatmap classification result
sc_hmClass(crParkall, sgrp, max=5000, isBig=TRUE, cCol=F, font=8)

stQuery <- get_cate(classRes = crParkall, sampTab = stQuery, dLevel = "seurat_clusters", sid = "sample_name", nrand = nqRand)

sc_violinClass(sampTab = stQuery, classRes = crParkall, sid = "sample_name", dLevel = "seurat_clusters", addRand = nqRand)

saveRDS(stQuery, 'stQuery_SingleCellNet_annotation_added.rds')



#------------------------------------------------------
# add the category column in the stQuery to object@meta.data
#------------------------------------------------------

rm(list = setdiff(ls(), "stQuery"))
object<- readRDS('object.final.rds')
metadata<- object@meta.data

identical(rownames(metadata), rownames(stQuery))

object@meta.data$category_SingleCellNet<- stQuery$category

saveRDS(object, 'object.final.rds')








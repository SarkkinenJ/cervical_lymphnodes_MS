library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)

# for annotation
library(SingleR)
library(RColorBrewer)
library(celldex)
library(scCATCH)

#### First analysis for each individual sample is done using analyzeSinglescRNAseqData.R from the shell with runSeuratInBatches.sh ####


#### Collect RDS files for the analyses of all samples, and combine the samples into one seurat object for analysis. ####

# The directory (/scratch/project_2005392/JoonaDawitAnalysisResults/) should be changed to the directory that holds the rds files of the analysis on the individual samples
seuratObjectsListNames = list.files("/scratch/project_2005392/JoonaDawitAnalysisResults/",pattern = "rds",full.names=T,recursive=T)


seuratObjectsList <- lapply(seuratObjectsListNames, function(x) readRDS(file = x))
names(seuratObjectsList) <- sapply(basename(seuratObjectsListNames),function(x) unlist(strsplit(x,"_"))[1])

patientStatus = c(rep("MS",5),rep("CTRL",3))
names(patientStatus) <- names(seuratObjectsList)

names(seuratObjectsList) <- sapply(1:length(patientStatus),function(x){ifelse(patientStatus[x]=="MS",names(patientStatus)[x],gsub("MS","CTRL",names(patientStatus)[x]))})


#### first simply merge the samples and evaluate the level of batch effect ####
msFNA.merged <- merge(seuratObjectsList[[1]], y = seuratObjectsList[-1], add.cell.ids = names(seuratObjectsList), project = "MSFNA")


# add the sample name and group status meta datas to the meta data of the merged object
msFNA.merged@meta.data$sampleName = sapply(rownames(msFNA.merged@meta.data),
                                           function(x) unlist(strsplit(x,"_"))[1])


msFNA.merged@meta.data$sampleGroup = "MS"
msFNA.merged@meta.data$sampleGroup[grepl("CTRL",msFNA.merged@meta.data$sampleName)] <- "CTRL"

# total number of cells form CTRL and MS samples
table(msFNA.merged@meta.data$sampleGroup)


# Run the standard workflow for visualization and clustering
msFNA.merged = NormalizeData(msFNA.merged, normalization.method = "LogNormalize", scale.factor = 10000)
msFNA.merged = FindVariableFeatures(msFNA.merged, selection.method = "vst", nfeatures = 3000)
msFNA.merged <- ScaleData(msFNA.merged, verbose = FALSE)
msFNA.merged <- RunPCA(msFNA.merged, verbose = FALSE) # variable features are used to run the pca by default


# pick PCs by JackStraw method
msFNA.merged <- JackStraw(msFNA.merged, prop.freq=0.1, dims = 50, num.replicate = 100)
msFNA.merged <- ScoreJackStraw(msFNA.merged, dims = 1:50)
JackStrawPlot(msFNA.merged, dims = 1:50)
ggsave(paste0("mergedData","_JackstrawPvalues.pdf"), height = 297, width = 210, units = "mm")


score.df <- JS(object = msFNA.merged[['pca']], slot = 'overall')
selectedPCs <- which(score.df[,2] < 0.05)



# t-SNE and Clustering
msFNA.merged <- RunUMAP(msFNA.merged, reduction = "pca", dims = selectedPCs)
msFNA.merged <- FindNeighbors(msFNA.merged, reduction = "pca", dims = selectedPCs)
msFNA.merged <- FindClusters(msFNA.merged, resolution = 1)



DimPlot(msFNA.merged, reduction = "umap", group.by = "sampleGroup")
ggsave(paste0("mergedData","clusteringByCondition.pdf"), height = 210, width = 320, units = "mm")

DimPlot(msFNA.merged, reduction = "umap", group.by = "sampleName")
ggsave(paste0("mergedData","clusteringBysampleName.pdf"), height = 210, width = 320, units = "mm")

DimPlot(msFNA.merged, reduction = "umap", label = TRUE, group.by=c("sampleGroup","SingleR.immune.monaco.mainLabels"), pt.size= 0.5)
ggsave(paste0("mergedData","clusteringByConditionAndAnnotation.pdf"), height = 210, width = 320, units = "mm")

DimPlot(msFNA.merged, reduction = "umap", label = TRUE, shape.by="sampleGroup",group.by=c("SingleR.immune.monaco.mainLabels"), pt.size= 0.5,label.size = 3)
ggsave(paste0("mergedData","clusteringSampleGroupAndAnnotation.pdf"), height = 210, width = 320, units = "mm")


# annotation after merging

# This function returns annotations form singleR, helps compare against the annotations done on the individual samples.
annotateCells <- function(srtObj){
  
  if(!exists("hpca.se"))
    hpca.se <- HumanPrimaryCellAtlasData()
  
  if(!exists("immune.monaco"))
    immune.monaco <- MonacoImmuneData()
  
  # annotate cells 
  
  normData <- GetAssayData(srtObj,assay="RNA",slot="data")
  
  #1
  hpcaMain <- SingleR(test = normData, ref = hpca.se, assay.type.test=1,
                      labels = hpca.se$label.main)$pruned.labels
  
  hpcaMain[is.na(hpcaMain)] <- "unannotated"
  
  #2
  hpcaFine <- SingleR(test = normData, ref = hpca.se, assay.type.test=1,
                      labels = hpca.se$label.fine)$pruned.labels
  hpcaFine[is.na(hpcaFine)] <- "unannotated"
  
  
  #3
  immMonacoMain <- SingleR(test = normData, ref = immune.monaco, assay.type.test=1,
                           labels = immune.monaco$label.main)$pruned.labels
  
  immMonacoMain[is.na(immMonacoMain)] <- "unannotated"
  
  
  #4
  immMonacoFine <- SingleR(test = normData, ref = immune.monaco, assay.type.test=1,
                           labels = immune.monaco$label.fine)$pruned.labels
  
  immMonacoFine[is.na(immMonacoFine)] <- "unannotated"
  
  return(list(hcmain=hpcaMain,hcfine=hpcaFine,immMain=immMonacoMain,immFine=immMonacoFine))
  
}


annotationResuls = annotateCells(msFNA.merged) # we pass the merged seurat object to the annotation function.

# add annotations form the merged data to the meta data

msFNA.merged[["hpcaMain"]] <- annotationResuls[[1]]
msFNA.merged[["hpcaFine"]] <- annotationResuls[[2]]
msFNA.merged[["immMain"]] <- annotationResuls[[3]]
msFNA.merged[["immFine"]] <- annotationResuls[[4]]




## save merged data
saveRDS(msFNA.merged, file = paste0("MS5_FNA","_merged2.rds"))


##### Integration of the samples instead of merging for analysis #### 


seuratObjectsList.forIntegration <- SplitObject(msFNA.merged, split.by = "sampleName")

seuratObjectsList.forIntegration <- lapply(X = seuratObjectsList.forIntegration, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 3000)
})

features <- SelectIntegrationFeatures(object.list = seuratObjectsList.forIntegration,nfeatures = 3000)
length(features)

# integration
selected.anchors <- FindIntegrationAnchors(object.list = seuratObjectsList.forIntegration, anchor.features = features,reduction = "rpca")
msFNA.integrated <- IntegrateData(anchorset = selected.anchors)


# Run the standard workflow for visualization and clustering
msFNA.integrated <- ScaleData(msFNA.integrated, verbose = FALSE)
msFNA.integrated <- RunPCA(msFNA.integrated, npcs = 50, verbose = FALSE)


msFNA.integrated <- RunUMAP(msFNA.integrated, reduction = "pca", dims = 1:50)
msFNA.integrated <- FindNeighbors(msFNA.integrated, reduction = "pca", dims = 1:50)
msFNA.integrated <- FindClusters(msFNA.integrated, resolution = 1)


# Visualization
head(msFNA.integrated@meta.data)

DimPlot(msFNA.integrated, reduction = "umap", label = TRUE, shape.by="sampleGroup",group.by=c("sampleName","SingleR.immune.monaco.mainLabels"), pt.size= 0.5,label.size = 3)
ggsave(paste0("integratedData","clusteringSampleGroupAndAnnotation.pdf"), height = 210, width = 320, units = "mm")


DimPlot(msFNA.integrated, reduction = "umap", label = TRUE, shape.by="sampleGroup",group.by=c("sampleName","SingleR.HumanPrimaryCellAtlasData.mainLabels"), pt.size= 0.5,label.size = 3)
ggsave(paste0("integratedData","clusteringSampleGroupAndHPCAAnnotation.pdf"), height = 210, width = 320, units = "mm")



DimPlot(msFNA.integrated, reduction = "umap", split.by = "sampleGroup")
ggsave(paste0("integratedData","clusteringBySampleGroup.pdf"), height = 210, width = 320, units = "mm")

msFNA.integrated@meta.data$SingleR.immune.monaco.mainLabels[is.na(msFNA.integrated@meta.data$SingleR.immune.monaco.mainLabels)] <- "unannotated"


DefaultAssay(msFNA.integrated) <- "RNA"


saveRDS(msFNA.integrated, file = paste0("MS5_FNA","_integrated2.rds"))









## Comment: we decided to use the merged version instead of the integrated version since the level of batch effect was minimal. Cells clustered by cell type
## rather than any other variable such as by sample or patient status.
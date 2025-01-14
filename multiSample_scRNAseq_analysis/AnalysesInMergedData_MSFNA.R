library(Seurat)
library(patchwork)
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(reshape2)

# for pseudobulking
library(muscat)

# for pathway analysis
library(pathifier)
library(fgsea)
library(msigdbr)
library(EnhancedVolcano)
require(graph)
require(ROntoTools) # topological pathway analysis with pathwayExpress
library(org.Hs.eg.db)

# for annotation
library(SingleR)
library(RColorBrewer)
library(celldex)
library(scCATCH)



# Read in merged data ####
msFNAmerged <- readRDS(file = "MS5_FNA_merged2.rds")

#### Evaluate batch effect visually ####
# There cell types appear to be well aligned between samples and sample groups in the merged data, with no significant batch effect

DimPlot(msFNAmerged, reduction = "umap", group.by = "sampleGroup")
DimPlot(msFNAmerged, reduction = "umap", group.by = "sampleName")
DimPlot(msFNAmerged, reduction = "umap", label = TRUE, group.by=c("sampleGroup","immMain"), pt.size= 0.5)
DimPlot(msFNAmerged, reduction = "umap", label = TRUE, group.by=c("sampleName","immMain"), pt.size= 0.5)

# total number of cells
totalNumberOfCells <- nrow(msFNAmerged@meta.data)

numberOfCellsPersample <- table(msFNAmerged@meta.data$sampleName)

par(mar = c(7, 4.1, 4.1, 2.1))
barplot(100 * (sort(table(msFNAmerged@meta.data$immMain),decreasing=T)/sum(table(msFNAmerged@meta.data$immMain))),
        las = 2,
        ylim=c(0,60),ylab="Detected cell types (%)")

ggsave("mergedMS5FNA_CellTypesProportion.pdf")

# number of detected cell types
write.csv(sort(table(msFNAmerged@meta.data$immMain),decreasing=T),file="mergedMS5FNA_CellTypesCounts.csv")
par(mar = c(7, 4.1, 4.1, 2.1))
barplot(sort(table(msFNAmerged@meta.data$immMain),decreasing=T),
        las = 2,
        ylim=c(0,50000),ylab="Detected cell types")

ggsave("mergedMS5FNA_CellTypesCounts.pdf")



# To check the batch effect further, we check the existence of batch effect within a single cluster of cells.

msFNAmerged.cd4subset <- subset(x = msFNAmerged, subset = immMain == "CD4+ T cells")

DimPlot(msFNAmerged.cd4subset, reduction = "umap", label = TRUE, shape.by="sampleGroup",group.by=c("sampleName","immMain"), pt.size= 0.5,label.size = 3)
ggsave(paste0("mergedData","clusteringSampleGroupAndAnnotation_CD4cellsOnly.pdf"), height = 210, width = 320, units = "mm")


# Run the standard workflow for visualization and clustering
msFNAmerged.cd4subset = NormalizeData(msFNAmerged.cd4subset, normalization.method = "LogNormalize", scale.factor = 10000)
msFNAmerged.cd4subset = FindVariableFeatures(msFNAmerged.cd4subset, selection.method = "vst", nfeatures = 3000)
msFNAmerged.cd4subset <- ScaleData(msFNAmerged.cd4subset, verbose = FALSE)
msFNAmerged.cd4subset <- RunPCA(msFNAmerged.cd4subset, verbose = FALSE) # variable features are used to run the pca by default


# t-SNE and Clustering
msFNAmerged.cd4subset <- RunUMAP(msFNAmerged.cd4subset, reduction = "pca", dims = 1:30)
msFNAmerged.cd4subset <- FindNeighbors(msFNAmerged.cd4subset, reduction = "pca", dims = 1:30)
msFNAmerged.cd4subset <- FindClusters(msFNAmerged.cd4subset, resolution = 1)


# Visualization
DimPlot(msFNAmerged.cd4subset, reduction = "umap", group.by = "sampleGroup")
ggsave(paste0("mergedData","clusteringSampleGroup_CD4cellsOnlyAfterReclustering.pdf"), height = 210, width = 320, units = "mm")


col_vector = unlist(mapply(brewer.pal, c(8,9,8,9), c("Accent","Set1","Dark2","Pastel1")))

DimPlot(msFNAmerged.cd4subset, reduction = "umap", group.by = "sampleName",label=T, cols=col_vector[1:length(unique(msFNAmerged.cd4subset@meta.data$sampleName))])
ggsave(paste0("combinedAnalysis/mergedData","clusteringSampleName_CD4cellsOnlyAfterReclustering.pdf"), height = 210, width = 320, units = "mm")

DimPlot(msFNAmerged.cd4subset, reduction = "umap", group.by = "sampleName",label=T, split.by="sampleGroup",cols=col_vector[1:length(unique(msFNAmerged.cd4subset@meta.data$sampleName))])
ggsave(paste0("combinedAnalysis/mergedData","clusteringSampleName2_CD4cellsOnlyAfterReclustering.pdf"), height = 210, width = 320, units = "mm")

DimPlot(msFNAmerged.cd4subset, reduction = "umap", label = TRUE, group.by=c("sampleGroup","immMain"), pt.size= 0.5)


DimPlot(msFNAmerged.cd4subset, reduction = "umap", label = TRUE, shape.by="sampleGroup",group.by=c("sampleName","immMain"), pt.size= 0.5,label.size = 3)
ggsave(paste0("mergedData","clusteringSampleGroupAndAnnotation_CD4cellsOnlyAfterReclustering.pdf"), height = 210, width = 320, units = "mm")

## It looks like there is not much batch effect even in cell subpopulations, e.g like in the CD4+ cells, no sample-specific or group-specific
## batch effect is seen.



#### Read in citeseq adt data for the merged dataset ####
samplesWithNoCiteSeq <- c()
AdtDataList <- list()

for(sm in unique(msFNAmerged@meta.data$orig.ident)){
  dataSample = paste0("/scratch/project_2005392/Joona/MS_FNA/scRNAseq/",sm)
  samp <- Read10X(dataSample)
  
  if(!is.list(samp)){
    #samp = samp$`Gene Expression`
    samplesWithNoCiteSeq <- c(samplesWithNoCiteSeq,sm)
  }else{
    citeseqCounts <- samp$`Antibody Capture`
    rownames(citeseqCounts) <- sapply(rownames(citeseqCounts),function(x) unlist(strsplit(x,"_"))[1])
    adt_assay <- CreateAssayObject(counts = citeseqCounts)
    AdtDataList[[sm]] <- adt_assay
    AdtDataList[[sm]] <- RenameCells(object = AdtDataList[[sm]], new.names = paste0(sm,"_",colnames(AdtDataList[[sm]])))
  }

}


msFNAmerged.samsWithCiteseq <- subset(x = msFNAmerged, subset = orig.ident %in% names(AdtDataList))


# combine all the adt data of all samples into one data matrix
combinedAdtData = NULL

for(sn in names(AdtDataList)){
  colidx <- which(colnames(AdtDataList[[sn]]) %in% rownames(msFNAmerged.samsWithCiteseq@meta.data[msFNAmerged.samsWithCiteseq@meta.data$orig.ident==sn,]))
  
  if(length(colidx) == 0){
    # this is a control sample. Change the cell prefix to control
    tempadt <- GetAssayData(object = AdtDataList[[sn]], slot = "data")
    colnames(tempadt) <- gsub("MS","CTRL",colnames(AdtDataList[[sn]]))
    
    colidx <- which(colnames(tempadt) %in% rownames(msFNAmerged.samsWithCiteseq@meta.data[msFNAmerged.samsWithCiteseq@meta.data$orig.ident==sn,]))
    adt_cellsToKeep <- colnames(tempadt)[colidx]
    tempadt <- tempadt[,adt_cellsToKeep]
    combinedAdtData <- cbind(combinedAdtData,tempadt)
  }else{
    adt_cellsToKeep <- colnames(AdtDataList[[sn]])[colidx]
    tempadt <- GetAssayData(object = AdtDataList[[sn]], slot = "data")
    tempadt <- tempadt[,adt_cellsToKeep]
    combinedAdtData <- cbind(combinedAdtData,tempadt)
  }
    print(sn)
    print(length(adt_cellsToKeep))
}


# Add the combined adt data as an assay to the merged seurat object
combined_adt_assay <- CreateAssayObject(counts = combinedAdtData)
msFNAmerged.samsWithCiteseq[["ADT"]] <- combined_adt_assay

DimPlot(msFNAmerged.samsWithCiteseq, label = TRUE,group.by=c("immMain","sampleGroup"))
DimPlot(msFNAmerged.samsWithCiteseq, label = TRUE,group.by=c("immMain","sampleName"))
ggsave("msFNAmerged.samsWithCiteseq_RNAseqclusters.pdf", height = 210, width = 320, units = "mm")



msFNAmerged.samsWithCiteseq <- NormalizeData(msFNAmerged.samsWithCiteseq, assay = "ADT", normalization.method = "CLR")
msFNAmerged.samsWithCiteseq <- ScaleData(msFNAmerged.samsWithCiteseq, assay = "ADT")

# do the clustering on the adt data

DefaultAssay(msFNAmerged.samsWithCiteseq) <- "ADT"
msFNAmerged.samsWithCiteseq <- RunPCA(msFNAmerged.samsWithCiteseq, features = rownames(msFNAmerged.samsWithCiteseq), reduction.name = "pca_adt", reduction.key = "pca_adt_", 
               verbose = FALSE)
# pca plot
DimPlot(msFNAmerged.samsWithCiteseq, reduction = "pca_adt",group.by=c("immMain","sampleGroup"))

ElbowPlot(msFNAmerged.samsWithCiteseq,ndims = 50)

# find the clusters in adt
msFNAmerged.samsWithCiteseq[["rnaClusterID"]] <- Idents(msFNAmerged.samsWithCiteseq) # save RNA clustering result

msFNAmerged.samsWithCiteseq <- FindNeighbors(msFNAmerged.samsWithCiteseq, dims = 1:30, reduction="pca_adt")
msFNAmerged.samsWithCiteseq <- FindClusters(msFNAmerged.samsWithCiteseq,graph.name="ADT_nn")
msFNAmerged.samsWithCiteseq.markers <- FindAllMarkers(msFNAmerged.samsWithCiteseq, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
msFNAmerged.samsWithCiteseq.markers.top <- msFNAmerged.samsWithCiteseq.markers %>%
  group_by(cluster) %>%
  slice_max(n = 2, order_by = avg_log2FC)

write.csv(msFNAmerged.samsWithCiteseq.markers,file="msFNAmerged.samsWithCiteseq.markers.csv")




#### Annotate cells with citeseq data and manual inspection ####

annotateCells <- function(srtObj){
  
  if(!exists("hpca.se"))
    hpca.se <- HumanPrimaryCellAtlasData()
  
  if(!exists("immune.monaco"))
    immune.monaco <- MonacoImmuneData()
  
  # annotate cells 
  
  normData <- GetAssayData(srtObj,assay="ADT",slot="data")
  
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


adtAnnotationResuls = annotateCells(msFNAmerged.samsWithCiteseq) # we pass the merged seurat object to the annotation function.

# add annotations form the merged data to the meta data
msFNAmerged.samsWithCiteseq[["adt_hpcaMain"]] <- adtAnnotationResuls[[1]]
msFNAmerged.samsWithCiteseq[["adt_hpcaFine"]] <- adtAnnotationResuls[[2]]
msFNAmerged.samsWithCiteseq[["adt_immMain"]] <- adtAnnotationResuls[[3]]
msFNAmerged.samsWithCiteseq[["adt_immFine"]] <- adtAnnotationResuls[[4]]


# ADT annotation of the unannotated population in the RNAseq data
msFNAmerged.samsWithCiteseq.unannotated <- msFNAmerged.samsWithCiteseq@meta.data[msFNAmerged.samsWithCiteseq@meta.data$immFine=="unannotated",]
adtAnnotationOfUnannotated <- t(table(msFNAmerged.samsWithCiteseq.unannotated$immFine,msFNAmerged.samsWithCiteseq.unannotated$adt_immFine))
write.csv(adtAnnotationOfUnannotated,file="adtAnnotationOfUnannotated.csv")


# adt clusters vs RNA cell annotation
adtclustervsRNAann <- table(msFNAmerged.samsWithCiteseq@meta.data$ADT_nn_res.0.8,msFNAmerged.samsWithCiteseq@meta.data$immFine)
adtclustersvsRNASingleR <- cbind(rownames(adtclustervsRNAann),apply(adtclustervsRNAann,1,function(x) colnames(adtclustervsRNAann)[which.max(x)]))

# adt clusters vs adt cell annotation
adtclustervsRNAann2 <- table(msFNAmerged.samsWithCiteseq@meta.data$ADT_nn_res.0.8,msFNAmerged.samsWithCiteseq@meta.data$adt_immFine)
adtclustersvsADTSingleR <- cbind(rownames(adtclustervsRNAann2),apply(adtclustervsRNAann2,1,function(x) colnames(adtclustervsRNAann2)[which.max(x)]))

adtclustersvsRNAandADTsingleRannotations <- cbind(adtclustersvsRNASingleR,adtclustersvsADTSingleR[,2])
colnames(adtclustersvsRNAandADTsingleRannotations) <- c("ADTcluster","RNAsingleRAnnotation","ADTsingleRAnnotation")
write.csv(adtclustersvsRNAandADTsingleRannotations,file="adtclustersvsRNAandADTsingleRannotations.csv")


## Annotating RNA clusters based on adt data

# RNA cluster markers, and annotations for each cluster

msFNAmerged.RNAdata.markers <- FindAllMarkers(msFNAmerged, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

# Find RNA markers for citeseq clusters
Idents(object = msFNAmerged.samsWithCiteseq) <- "ADT_nn_res.0.8"
DefaultAssay(msFNAmerged.samsWithCiteseq) <- "RNA"

msFNAmerged.samsWithCiteseq.citeseqClusterRNAmarkers <- FindAllMarkers(msFNAmerged.samsWithCiteseq, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)


msFNAmerged.samsWithCiteseq.unannotated <- msFNAmerged.samsWithCiteseq@meta.data[msFNAmerged.samsWithCiteseq@meta.data$immFine=="unannotated",]

# ADT annotation of the unannotated population in the RNAseq data
adtAnnotationOfUnannotated <- t(table(msFNAmerged.samsWithCiteseq.unannotated$immFine,msFNAmerged.samsWithCiteseq.unannotated$adt_immFine))
write.csv(adtAnnotationOfUnannotated,file="adtAnnotationOfUnannotated.csv")

writeCellTypePredictionsPerCluster <- function(seuratObj,clusterLabels="RNA_snn_res.1",annotationType="immFine",prefix="msFNAmerged.RNAseq"){
    for(unqType in unique(seuratObj[[clusterLabels]][,1])){
      #cellsToAdd = rownames(seuratObj@meta.data[msFNAmerged[[clusterLabels]][,1]==unqType,])
      tempData <- seuratObj@meta.data[seuratObj[[clusterLabels]][,1]==unqType,]
      AnnotationOfLabel <- t(table(tempData[,clusterLabels],tempData[,annotationType]))
      AnnotationOfLabel <- AnnotationOfLabel[,colnames(AnnotationOfLabel)==unqType,drop=F]
      
      AnnotationOfLabel <- AnnotationOfLabel[order(AnnotationOfLabel[,1],decreasing=T),,drop=F]
      
      #if(grepl("/",unqType)) unqType = gsub("/","-",unqType)
      
      write.csv(AnnotationOfLabel,file=paste0("/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/mergedClusterMarkers/",prefix,"_",unqType,"_",annotationType,"_annotationRank.csv"))
      }
  

}

writeCellTypePredictionsPerCluster(msFNAmerged,clusterLabels="RNA_snn_res.1",annotationType="immFine",prefix="msFNAmerged.RNAseq")

writeCellTypePredictionsPerCluster(msFNAmerged.samsWithCiteseq,clusterLabels="ADT_nn_res.0.8",annotationType="adt_immFine",prefix="msFNAmerged.Citeseq")

writeCellTypePredictionsPerCluster(msFNAmerged.samsWithCiteseq,clusterLabels="ADT_nn_res.0.8",annotationType="immFine",prefix="msFNAmerged.Citeseq")


DimPlot(msFNAmerged.samsWithCiteseq, label = TRUE, reduction = "umap_adt",group.by="ADT_nn_res.0.8")



# We use max adt annotation per RNA cluster and assign that as the cluster identity

rnaclustervsAdtAnn <- table(msFNAmerged.samsWithCiteseq@meta.data$RNA_snn_res.1,msFNAmerged.samsWithCiteseq@meta.data$adt_immFine)
rnaclustervsAdtRNASingleR <- cbind(rownames(rnaclustervsAdtAnn),apply(rnaclustervsAdtAnn,1,function(x) colnames(rnaclustervsAdtAnn)[which.max(x)]))
colnames(rnaclustervsAdtRNASingleR) <- c("rnaclusterid","maxAdtAnnotation")
rnaClusterAnnotationFromAdt <- rnaclustervsAdtRNASingleR[,2][match(as.character(msFNAmerged.samsWithCiteseq@meta.data$RNA_snn_res.1),rnaclustervsAdtRNASingleR[,1])]

rnaClusterAnnotationFromAdtAllData <- rnaclustervsAdtRNASingleR[,2][match(as.character(msFNAmerged@meta.data$RNA_snn_res.1),rnaclustervsAdtRNASingleR[,1])]


msFNAmerged.samsWithCiteseq[["rnaClusterAnnotationFromAdt"]] <- rnaClusterAnnotationFromAdt

msFNAmerged[["rnaClusterAnnotationFromAdt"]] <- rnaClusterAnnotationFromAdtAllData


## 

msFNAmerged.samsWithCiteseq <- RunUMAP(msFNAmerged.samsWithCiteseq, dims = 1:30,reduction = "pca_adt",reduction.name = "umap_adt",reduction.key = "UMAPadt_")

DimPlot(msFNAmerged.samsWithCiteseq, label = TRUE, reduction = "umap",group.by=c("immMain","adt_immMain"))
ggsave("msFNAmerged.samsWithCiteseq_RNAseqClustering_immMain_adt_annotation.pdf", height = 210, width = 320, units = "mm")

p1 <- DimPlot(msFNAmerged.samsWithCiteseq, label = T,repel=T,label.size=2.5, reduction = "umap",group.by="immFine") + NoLegend()
p2 <- DimPlot(msFNAmerged.samsWithCiteseq, label = T,repel=T,label.size=2,reduction = "umap",group.by="adt_immFine") + NoLegend()

p1|p2
ggsave("msFNAmerged.samsWithCiteseq_RNAseqClustering_immFine_adt_annotation.pdf", height = 210, width = 320, units = "mm")

# use the RNA cluster level annotation using the maximum adt_immFine values for each RNA cluster
p1 <- DimPlot(msFNAmerged.samsWithCiteseq, label = T,repel=T,label.size=2.5, reduction = "umap",group.by="RNA_snn_res.1") + NoLegend()
p2 <- DimPlot(msFNAmerged.samsWithCiteseq, label = T,repel=T,label.size=2,reduction = "umap",group.by="rnaClusterAnnotationFromAdt") + NoLegend()

p1|p2
ggsave("msFNAmerged.samsWithCiteseq_RNAseqClusters_AnnotatedByAdtannotation.pdf", height = 210, width = 320, units = "mm")

#
p1 <- DimPlot(msFNAmerged.samsWithCiteseq, label = T,repel=T,label.size=2.5, reduction = "umap",group.by="immFine") + NoLegend()
p2 <- DimPlot(msFNAmerged.samsWithCiteseq, label = T,repel=T,label.size=2,reduction = "umap",group.by="rnaClusterAnnotationFromAdt") + NoLegend()

p1|p2
ggsave("msFNAmerged.samsWithCiteseq_RNAseqClusters_AnnotatedByAdtannotation2.pdf", height = 210, width = 320, units = "mm")


DimPlot(msFNAmerged.samsWithCiteseq, label = TRUE, reduction = "umap_adt",group.by="adt_immMain")

DimPlot(msFNAmerged.samsWithCiteseq, label = TRUE,group.by=c("immMain","sampleGroup"))
DimPlot(msFNAmerged.samsWithCiteseq, label = TRUE,group.by=c("immMain","sampleName"))
DimPlot(msFNAmerged.samsWithCiteseq, label = TRUE,group.by=c("sampleGroup","sampleName"))
ggsave("msFNAmerged.samsWithCiteseq_citeseqClustering_umapPlot.pdf", height = 210, width = 320, units = "mm")
ggsave("msFNAmerged.samsWithCiteseq_citeseqClustering_umapPlot4.pdf", height = 210, width = 320, units = "mm")


DefaultAssay(msFNAmerged.samsWithCiteseq) <- "RNA"
DimPlot(msFNAmerged.samsWithCiteseq,  reduction = "umap",label = TRUE,group.by=c("immMain","sampleGroup"))



## Annotation of cells after manual inspection of the cluster markers detected from the RNA and ADT datasets

# Manual aided annotation of cells. We had a closer look at the adt clusters, and suggested most were good and some changes needed
# we followed manual guidance to re-annotate the cells (This was mostly done interactively)
# and assigned citeseq guided annotation to the citeseq clusters manually.

adtclustersvsAdtAnn <- table(msFNAmerged.samsWithCiteseq@meta.data$ADT_nn_res.0.8,msFNAmerged.samsWithCiteseq@meta.data$adt_immFine)

adtclustersvsAdtAnnMaxAnn <- cbind(rownames(adtclustersvsAdtAnn),apply(adtclustersvsAdtAnn,1,function(x) colnames(adtclustersvsAdtAnn)[which.max(x)]))
colnames(adtclustersvsAdtAnnMaxAnn) <- c("adtclustersid","maxAdtAnnotation")

# Joona's manual annotation.
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==0,2] <- "NaiveCD4.1"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==1,2] <- "NaiveCD4.2"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==2,2] <- "NaiveCD4.3"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==3,2] <- "Treg1"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==4,2] <- "Tfh.1"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==5,2] <- "NaiveCD8.1"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==6,2] <- "B cells"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==7,2] <- "NaiveCD4.4"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==8,2] <- "B cells"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==9,2] <- "Tfr"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==10,2] <- "NaiveCD4.5"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==11,2] <- "memCD8"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==12,2] <- "Tfh.2"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==13,2] <- "B cells"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==14,2] <- "NaiveCD8.2"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==15,2] <- "Myeloid"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==16,2] <- "Platelets"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==17,2] <- "NK cells"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==18,2] <- "NaiveCD4.6"
adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==19,2] <- "gdTcell"

# we merged the small adt cluster 20, to its closest adt cluster 5 and as both have been annotated as cd8, we merge them.

adtclustersvsAdtAnnMaxAnn[rownames(adtclustersvsAdtAnnMaxAnn)==20,2] <- "NaiveCD8.1"

# assign citeseq and manual annotation to data
adtClusterAnnotationFromAdtManual <- adtclustersvsAdtAnnMaxAnn[,2][match(as.character(msFNAmerged.samsWithCiteseq@meta.data$ADT_nn_res.0.8),adtclustersvsAdtAnnMaxAnn[,1])]
msFNAmerged.samsWithCiteseq[["adtClusterAnnotationFromAdtManual"]] <- adtClusterAnnotationFromAdtManual
DimPlot(msFNAmerged.samsWithCiteseq, label = TRUE, reduction = "umap",group.by=c("RNA_snn_res.1","adtClusterAnnotationFromAdtManual"))


# evaluate rna clusters and the adt manual labels
rnaclustersvsAdtManAnn <- table(msFNAmerged.samsWithCiteseq@meta.data$RNA_snn_res.1,msFNAmerged.samsWithCiteseq@meta.data$adtClusterAnnotationFromAdtManual)
rnaclustersvsAdtCl <- table(msFNAmerged.samsWithCiteseq@meta.data$RNA_snn_res.1,msFNAmerged.samsWithCiteseq@meta.data$ADT_nn_res.0.8)


# set all naivecd4.1 etc to naice CD4, and same for NaiveCD8s
adtClusterAnnotationFromAdtManual2 <- adtClusterAnnotationFromAdtManual

adtClusterAnnotationFromAdtManual2[grepl("NaiveCD4",adtClusterAnnotationFromAdtManual2)] <- "NaiveCD4"
adtClusterAnnotationFromAdtManual2[grepl("NaiveCD8",adtClusterAnnotationFromAdtManual2)] <- "NaiveCD8"

msFNAmerged.samsWithCiteseq[["adtClusterAnnotationFromAdtManual2"]] <- adtClusterAnnotationFromAdtManual2
DimPlot(msFNAmerged.samsWithCiteseq, label = TRUE, reduction = "umap",group.by=c("RNA_snn_res.1","adtClusterAnnotationFromAdtManual2"))


#Next assign the most frequent adt manual lables to RNA clusters (cluster level annotation)
rnaclustersvsAdtManAnn2 <- table(msFNAmerged.samsWithCiteseq@meta.data$RNA_snn_res.1,msFNAmerged.samsWithCiteseq@meta.data$adtClusterAnnotationFromAdtManual2)
rnaclustersvsAdtManAnn2MaxAnn <- cbind(rownames(rnaclustersvsAdtManAnn2),apply(rnaclustersvsAdtManAnn2,1,function(x) colnames(rnaclustersvsAdtManAnn2)[which.max(x)]))

# since some RNA clusters are appear very distinct in clusters but are given the same labels, we compare their markers and
# check if their annotation should be changed. 

# compare RNA cluster 15 from other naive CD4s, why is it so distinct ?
DefaultAssay(msFNAmerged.samsWithCiteseq) <- "RNA"
RNAcl.15.markers <- FindMarkers(msFNAmerged.samsWithCiteseq, 
                                ident.1 = "15", 
                                ident.2 = c("0","10","11"),
                                only.pos = T,
                                group.by="RNA_snn_res.1")

DefaultAssay(msFNAmerged.samsWithCiteseq) <- "ADT"
RNAcl.15.adtmarkers <- FindMarkers(msFNAmerged.samsWithCiteseq, logfc.threshold=0.1,
                                ident.1 = "15", 
                                ident.2 = c("0","10","11"),
                                only.pos = T,
                                group.by="RNA_snn_res.1")


RNAcl.15.markers2 <- RNAcl.15.markers[RNAcl.15.markers$p_val_adj < 0.05,] 
RNAcl.15.markers2 <- RNAcl.15.markers2[order(RNAcl.15.markers2$avg_log2FC,decreasing=T),]

write.csv(RNAcl.15.markers2,file="RNAcl.15vs0.10.11.markers2.csv")
write.csv(RNAcl.15.adtmarkers,file="RNAcl.15vs0.10.11.adtmarkers.csv")


# A main RNA markers in cluster 15 is KLRB1 and it is in the variable genes used for the RNA clustering.

"KLRB1" %in% VariableFeatures(msFNAmerged.samsWithCiteseq)
"KLRG1" %in% VariableFeatures(msFNAmerged.samsWithCiteseq)


DefaultAssay(msFNAmerged.samsWithCiteseq) <- "ADT"
Idents(msFNAmerged.samsWithCiteseq) <- "ADT_nn_res.0.8"

FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("KLRB1"))

FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("CD8A"),label=T,split.by = "sampleGroup")

Idents(msFNAmerged.samsWithCiteseq) <- "RNA_snn_res.1"
FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("KLRB1"),label=T,split.by = "sampleGroup")


# cluster 13 

DefaultAssay(msFNAmerged.samsWithCiteseq) <- "RNA"
RNAcl.13.markers <- FindMarkers(msFNAmerged.samsWithCiteseq, 
                                ident.1 = "13", 
                                ident.2 = c("0","10","11"),
                                only.pos = T,
                                group.by="RNA_snn_res.1")

DefaultAssay(msFNAmerged.samsWithCiteseq) <- "ADT"
RNAcl.13.adtmarkers <- FindMarkers(msFNAmerged.samsWithCiteseq, logfc.threshold=0.1,
                                   ident.1 = "13", 
                                   ident.2 = c("0","10","11"),
                                   only.pos = T,
                                   group.by="RNA_snn_res.1")


RNAcl.13.markers2 <- RNAcl.13.markers[RNAcl.13.markers$p_val_adj < 0.05,] 
RNAcl.13.markers2 <- RNAcl.13.markers2[order(RNAcl.13.markers2$avg_log2FC,decreasing=T),]

write.csv(RNAcl.13.markers2,file="RNAcl.13vs0.10.11.markers2.csv")
write.csv(RNAcl.13.adtmarkers,file="RNAcl.13vs0.10.11.adtmarkers.csv")


Idents(msFNAmerged.samsWithCiteseq) <- "RNA_snn_res.1"
FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("XIST"),label=T)

FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("NEAT1"),label=T)

#FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("IL7R"),label=T)



# differences between RNA clusters 7(predicted treg, but confusing) and 8(treg), and 12 (tfr)


FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("FOXP3"),label=T,split.by = "sampleGroup")

# so 8 is treg, and 12 is tfr. 7 is then more likely to be tfh, but what makes it different from the other tfhs like 1

## 7 vs 1 (another tfh)

DefaultAssay(msFNAmerged.samsWithCiteseq) <- "RNA"
RNAcl.7vs1.markers <- FindMarkers(msFNAmerged.samsWithCiteseq, 
                                ident.1 = "7", 
                                ident.2 = "1",
                                only.pos = T,
                                group.by="RNA_snn_res.1")


RNAcl.7vs1.markers2 <- RNAcl.7vs1.markers[RNAcl.7vs1.markers$p_val_adj < 0.05,] 
RNAcl.7vs1.markers2 <- RNAcl.7vs1.markers2[order(RNAcl.7vs1.markers2$avg_log2FC,decreasing=T),]

# 7 vs 0 (naive cd4)
RNAcl.7vs0.markers <- FindMarkers(msFNAmerged.samsWithCiteseq, 
                                  ident.1 = "7", 
                                  ident.2 = "0",
                                  only.pos = T,
                                  group.by="RNA_snn_res.1")


RNAcl.7vs0.markers2 <- RNAcl.7vs0.markers[RNAcl.7vs0.markers$p_val_adj < 0.05,] 
RNAcl.7vs0.markers2 <- RNAcl.7vs0.markers2[order(RNAcl.7vs0.markers2$avg_log2FC,decreasing=T),]


FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("KLF2"),label=T,split.by = "sampleGroup")

"CXCR5" %in% rownames(RNAcl.7vs1.markers2)
"PDCD1" %in% rownames(RNAcl.7vs1.markers2)

FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("CXCR5"),label=T,split.by = "sampleGroup")

FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("AHNAK"),label=T,split.by = "sampleGroup")


## For the b cell subtypes:
# cluster 17 is plasmablasts, 21 and 23 could be GC.4 and 5 are naive bcells, what is 26 & 30 ?
DefaultAssay(msFNAmerged.samsWithCiteseq) <- "RNA"
RNAcl.21vs4.5.markers <- FindMarkers(msFNAmerged.samsWithCiteseq, 
                                  ident.1 = "21", 
                                  ident.2 = c("4","5"),
                                  only.pos = T,
                                  group.by="RNA_snn_res.1")


RNAcl.21vs4.5.markers2 <- RNAcl.21vs4.5.markers[RNAcl.21vs4.5.markers$p_val_adj < 0.05,] 
RNAcl.21vs4.5.markers2 <- RNAcl.21vs4.5.markers2[order(RNAcl.21vs4.5.markers2$avg_log2FC,decreasing=T),]

FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("HMGB2"),label=T,split.by = "sampleGroup")
FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("BCL6"),label=T)
FeaturePlot(msFNAmerged.samsWithCiteseq, features = c("AICDA"),label=T)



"BCL6" %in% rownames(RNAcl.21vs4.5.markers2)
"AICDA" %in% rownames(RNAcl.21vs4.5.markers2)

RNAcl.21vs4.5.markers2[c("BCL6","AICDA"),]

# reclustering of the b cell populations only, only the bcells: cluster 5 (naive b cells), 4 (memory b cells)
# 5,4, 26,30,21,23,17
# We more or less know 23, no need to reclustering.
# determine 26, 30 ,very small populations, we can call them bcells for now (memory or naive).

# For B cells, 4 is memory, 5 is naive, 21 and 23 are GC,26 and 30 are small just check if they belong to 4 or 5 and assign them to either
# determine 13 and 15, 7 is tfh (not treg) and to be named tfh3, 1 is tfh1 and 22 is tfh2.


# Change the labels, make adjustments

rnaclustersvsAdtManAnn2MaxAnn[,2][rnaclustersvsAdtManAnn2MaxAnn[,1]==7] <- "Tfh.3"
rnaclustersvsAdtManAnn2MaxAnn[,2][rnaclustersvsAdtManAnn2MaxAnn[,1]==22] <- "Tfh.2"

rnaclustersvsAdtManAnn2MaxAnn[,2][rnaclustersvsAdtManAnn2MaxAnn[,1]==4] <- "Mem B cells"
rnaclustersvsAdtManAnn2MaxAnn[,2][rnaclustersvsAdtManAnn2MaxAnn[,1]==5] <- "Naive B cells"
rnaclustersvsAdtManAnn2MaxAnn[,2][rnaclustersvsAdtManAnn2MaxAnn[,1]==26] <- "B cells"
rnaclustersvsAdtManAnn2MaxAnn[,2][rnaclustersvsAdtManAnn2MaxAnn[,1]==30] <- "B cells"

rnaclustersvsAdtManAnn2MaxAnn[,2][rnaclustersvsAdtManAnn2MaxAnn[,1]==21] <- "GC B cells"
rnaclustersvsAdtManAnn2MaxAnn[,2][rnaclustersvsAdtManAnn2MaxAnn[,1]==23] <- "GC B cells"

rnaclustersvsAdtManAnn2MaxAnn[,2][rnaclustersvsAdtManAnn2MaxAnn[,1]==17] <- "Plasmablasts"

rnaclustersvsAdtManAnn2MaxAnn2 <- rnaclustersvsAdtManAnn2MaxAnn


# We can use the rnaclustersvsAdtManAnn2MaxAnn2 labels as a final label of the clusters.


# Next B cells: first recluster b cells, and redefine the clusters.

bcellClusters <- c("17",rnaclustersvsAdtManAnn2MaxAnn2[,1][grepl("B cells",rnaclustersvsAdtManAnn2MaxAnn2[,2])])

msFNAmerged.samsWithCiteseq.bcells <- subset(x = msFNAmerged.samsWithCiteseq, subset = RNA_snn_res.1 %in% bcellClusters)
msFNAmerged.samsWithCiteseq.bcells <- NormalizeData(msFNAmerged.samsWithCiteseq.bcells)
msFNAmerged.samsWithCiteseq.bcells <- FindVariableFeatures(msFNAmerged.samsWithCiteseq.bcells, selection.method = "vst", nfeatures = 3000)

length(intersect(VariableFeatures(msFNAmerged.samsWithCiteseq.bcells),VariableFeatures(msFNAmerged.samsWithCiteseq)))
length(setdiff(VariableFeatures(msFNAmerged.samsWithCiteseq.bcells),VariableFeatures(msFNAmerged.samsWithCiteseq)))

variableGenesInBcellsAndTotal = cbind(head(VariableFeatures(msFNAmerged.samsWithCiteseq.bcells),30),head(VariableFeatures(msFNAmerged.samsWithCiteseq),30))

all.genes.bcells <- rownames(msFNAmerged.samsWithCiteseq.bcells)
msFNAmerged.samsWithCiteseq.bcells <- ScaleData(msFNAmerged.samsWithCiteseq.bcells, features = all.genes.bcells)

s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
msFNAmerged.samsWithCiteseq.bcells <- CellCycleScoring(msFNAmerged.samsWithCiteseq.bcells, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)


msFNAmerged.samsWithCiteseq.bcells <- RunPCA(msFNAmerged.samsWithCiteseq.bcells, features = VariableFeatures(object = msFNAmerged.samsWithCiteseq.bcells))


msFNAmerged.samsWithCiteseq.bcells <- FindNeighbors(msFNAmerged.samsWithCiteseq.bcells, dims = 1:15)
msFNAmerged.samsWithCiteseq.bcells <- FindClusters(msFNAmerged.samsWithCiteseq.bcells, resolution = 0.1)
msFNAmerged.samsWithCiteseq.bcells <- RunUMAP(msFNAmerged.samsWithCiteseq.bcells, dims = 1:15)

DimPlot(msFNAmerged.samsWithCiteseq.bcells, reduction = "umap",label=T)
ggsave("/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/mergedBcellReclusteringMarkers/UmapPlot_msFNAmerged.samsWithCiteseq.bcells.pdf")

rnaCluster_bcellRNAreclusters <- table(msFNAmerged.samsWithCiteseq.bcells@meta.data$RNA_snn_res.1,msFNAmerged.samsWithCiteseq.bcells@meta.data$RNA_snn_res.0.1)
rnaCluster_bcellRNAreclusters <- rnaCluster_bcellRNAreclusters[rownames(rnaCluster_bcellRNAreclusters) %in% bcellClusters,]

write.csv(rnaCluster_bcellRNAreclusters,"/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/mergedBcellReclusteringMarkers/rnaCluster_bcellRNAreclusters.csv")

bcellRNAreclusters_immFine <- t(table(msFNAmerged.samsWithCiteseq.bcells@meta.data$RNA_snn_res.0.1,msFNAmerged.samsWithCiteseq.bcells@meta.data$immFine))
write.csv(bcellRNAreclusters_immFine,"bcellRNAreclusters_immFine.csv")

msFNAmerged.samsWithCiteseq.bcells.markers <- FindAllMarkers(msFNAmerged.samsWithCiteseq.bcells, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(msFNAmerged.samsWithCiteseq.bcells.markers,file="/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/mergedBcellReclusteringMarkers/mergedBcellReclutseringMarkers.csv")

jBregIntGenes <- c("IL10", "Lag3", "TIM1", "AHR", "MAPK3")
#
inDE.bmarkers <- jBregIntGenes %in% msFNAmerged.samsWithCiteseq.bcells.markers$gene

msFNAmerged.samsWithCiteseq.bcells.markers[msFNAmerged.samsWithCiteseq.bcells.markers$gene == "AHR", ]

FeaturePlot(msFNAmerged.samsWithCiteseq.bcells, features = c("BCL6"),label=T)
FeaturePlot(msFNAmerged.samsWithCiteseq.bcells, features = c("AHR"),label=T)
FeaturePlot(msFNAmerged.samsWithCiteseq.bcells, features = c("MAPK3"),label=T)
FeaturePlot(msFNAmerged.samsWithCiteseq.bcells, features = c("AICDA","BCL6","BCL2A1"),label=T)
ggsave("/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/mergedBcellReclusteringMarkers/FeaturePlotBCL6_BCL2A1_AICDA_msFNAmerged.samsWithCiteseq.bcells.pdf")

FeaturePlot(msFNAmerged.samsWithCiteseq.bcells, features = c("IL10"),label=T)

jBregIntGenes[jBregIntGenes %in% VariableFeatures(msFNAmerged.samsWithCiteseq.bcells)]


VlnPlot(msFNAmerged.samsWithCiteseq.bcells, features = c("S.Score"))
VlnPlot(msFNAmerged.samsWithCiteseq.bcells, features = c("G2M.Score"))
ggsave("/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/mergedBcellReclusteringMarkers/VlnPlot_G2M.score_msFNAmerged.samsWithCiteseq.bcells.pdf")


# recluster 3
msFNAmerged.samsWithCiteseq.bcells3 <- subset(x = msFNAmerged.samsWithCiteseq.bcells, subset = RNA_snn_res.0.1 == "3")

DefaultAssay(msFNAmerged.samsWithCiteseq.bcells3)

msFNAmerged.samsWithCiteseq.bcells3 <- NormalizeData(msFNAmerged.samsWithCiteseq.bcells3)
msFNAmerged.samsWithCiteseq.bcells3 <- FindVariableFeatures(msFNAmerged.samsWithCiteseq.bcells3, selection.method = "vst", nfeatures = 2000)

length(intersect(VariableFeatures(msFNAmerged.samsWithCiteseq.bcells),VariableFeatures(msFNAmerged.samsWithCiteseq.bcells3)))
length(setdiff(VariableFeatures(msFNAmerged.samsWithCiteseq.bcells3),VariableFeatures(msFNAmerged.samsWithCiteseq.bcells)))

all.genes.bcells3 <- rownames(msFNAmerged.samsWithCiteseq.bcells3)
msFNAmerged.samsWithCiteseq.bcells3 <- ScaleData(msFNAmerged.samsWithCiteseq.bcells3, features = all.genes.bcells3)

msFNAmerged.samsWithCiteseq.bcells3 <- RunPCA(msFNAmerged.samsWithCiteseq.bcells3, features = VariableFeatures(object = msFNAmerged.samsWithCiteseq.bcells3))

msFNAmerged.samsWithCiteseq.bcells3 <- FindNeighbors(msFNAmerged.samsWithCiteseq.bcells3, dims = 1:20)
msFNAmerged.samsWithCiteseq.bcells3 <- FindClusters(msFNAmerged.samsWithCiteseq.bcells3, resolution = 0.4)
msFNAmerged.samsWithCiteseq.bcells3 <- RunUMAP(msFNAmerged.samsWithCiteseq.bcells3, dims = 1:20)

DimPlot(msFNAmerged.samsWithCiteseq.bcells3, reduction = "umap",label=T)
ggsave("/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/mergedBcellReclusteringMarkers/UmapPlot_msFNAmerged.samsWithCiteseq.bcells.cluster3.pdf")

msFNAmerged.samsWithCiteseq.bcells3.markers <- FindAllMarkers(msFNAmerged.samsWithCiteseq.bcells3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(msFNAmerged.samsWithCiteseq.bcells3.markers,file="/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/mergedBcellReclusteringMarkers/msFNAmerged.samsWithCiteseq.bcells_cluster3_subclusterMarkers.csv")

VlnPlot(msFNAmerged.samsWithCiteseq.bcells3, features = c("S.Score"))
VlnPlot(msFNAmerged.samsWithCiteseq.bcells3, features = c("G2M.Score"))
ggsave("/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/mergedBcellReclusteringMarkers/VlnPlot_S.score_msFNAmerged.samsWithCiteseq.bcells3_subclusters.pdf")



FeaturePlot(msFNAmerged.samsWithCiteseq.bcells3, features = c("AICDA","BCL6","BCL2A1"),label=T)


# first make minor annotation changes.
rnaclustersvsAdtManAnn2MaxAnn2[,2][rnaclustersvsAdtManAnn2MaxAnn2[,1]==8] <- "Treg"
rnaclustersvsAdtManAnn2MaxAnn2[,2][rnaclustersvsAdtManAnn2MaxAnn2[,1]==30] <- "Mem B cells"
rnaclustersvsAdtManAnn2MaxAnn2[,2][rnaclustersvsAdtManAnn2MaxAnn2[,1]==26] <- "Mem B cells"


# then add the final annotation to the metadata of all samples and citeseq having samples, and save the workspaces

colnames(rnaclustersvsAdtManAnn2MaxAnn2) <- c("rnaclusterid","rnaclustersvsAdtManAnn2MaxAnn2")
rnaclustersvsAdtManAnn2MaxAnn2PerCell <- rnaclustersvsAdtManAnn2MaxAnn2[,2][match(as.character(msFNAmerged.samsWithCiteseq@meta.data$RNA_snn_res.1),rnaclustersvsAdtManAnn2MaxAnn2[,1])]

rnaclustersvsAdtManAnn2MaxAnn2AllData <- rnaclustersvsAdtManAnn2MaxAnn2[,2][match(as.character(msFNAmerged@meta.data$RNA_snn_res.1),rnaclustersvsAdtManAnn2MaxAnn2[,1])]

msFNAmerged.samsWithCiteseq[["rnaclustersvsAdtManAnn2MaxAnn2"]] <- rnaclustersvsAdtManAnn2MaxAnn2PerCell
msFNAmerged[["rnaclustersvsAdtManAnn2MaxAnn2"]] <- rnaclustersvsAdtManAnn2MaxAnn2AllData


# write the seurat objects

# keep samples with good number of cells, and remove samples with low number of cells
samnames = unique(msFNAmerged@meta.data$sampleName)
samnamesSelected = samnames[!samnames %in% c("MS001","CTRL012","CTRL013","CTRL014")]
msFNAmerged.goodSamples <- subset(x = msFNAmerged, subset = sampleName %in% samnamesSelected)


#### read in Final annotation table and annotate the clusters ####

# Final updated annotation of the cells using manual annotation (after several rounds of subclustering and interactive inspection of clusters and their markers from RNA and ADT data)

RNAcluster_finalAnnotation = read.table(file="RNAcluster_FinalAnnotation.txt",header=T,sep="\t")

# for the entire merged dataset: msFNAmerged
colnames(msFNAmerged@meta.data)[colnames(msFNAmerged@meta.data)=="rnaclustersvsAdtManAnn2MaxAnn2"] <- "rnaclustersvsAdtManAnn2MaxAnn2_old"
msFNAmerged[["rnaclustersvsAdtManAnn2MaxAnn2"]] <- RNAcluster_finalAnnotation$Final_Annotation[match(msFNAmerged@meta.data$seurat_clusters,RNAcluster_finalAnnotation$RNAcluster)]

# for the merged dataset that are good sized samples which are finally used for analysis and paper: msFNAmerged.goodSamples
colnames(msFNAmerged.goodSamples@meta.data)[colnames(msFNAmerged.goodSamples@meta.data)=="rnaclustersvsAdtManAnn2MaxAnn2"] <- "rnaclustersvsAdtManAnn2MaxAnn2_old"
msFNAmerged.goodSamples[["rnaclustersvsAdtManAnn2MaxAnn2"]] <- RNAcluster_finalAnnotation$Final_Annotation[match(msFNAmerged.goodSamples@meta.data$seurat_clusters,RNAcluster_finalAnnotation$RNAcluster)]

# for the merged dataset that has samples with citeseq only: msFNAmerged.samsWithCiteseq
# Also rename PTPRC, PTPRC.1 and PTPRC.2 gene names in the citeseq data to CD45RA,CD45RO & CD45 respectively
colnames(msFNAmerged.samsWithCiteseq@meta.data)[colnames(msFNAmerged.samsWithCiteseq@meta.data)=="rnaclustersvsAdtManAnn2MaxAnn2"] <- "rnaclustersvsAdtManAnn2MaxAnn2_old"
msFNAmerged.samsWithCiteseq[["rnaclustersvsAdtManAnn2MaxAnn2"]] <- RNAcluster_finalAnnotation$Final_Annotation[match(msFNAmerged.samsWithCiteseq@meta.data$RNA_snn_res.1,RNAcluster_finalAnnotation$RNAcluster)]

renameGenesInSeuratADT <- function(seuratObject,originalName,newName){
  
  rownames(seuratObject@assays[["ADT"]]@counts)[rownames(seuratObject@assays[["ADT"]]@counts) == originalName] <- newName
  rownames(seuratObject@assays[["ADT"]]@data)[rownames(seuratObject@assays[["ADT"]]@data) == originalName] <- newName
  rownames(seuratObject@assays[["ADT"]]@scale.data)[rownames(seuratObject@assays[["ADT"]]@scale.data) == originalName] <- newName
  
  rownames(seuratObject@assays[["ADT"]]@meta.features)[rownames(seuratObject@assays[["ADT"]]@meta.features) == originalName] <- newName
  
  return(seuratObject)
}

msFNAmerged.samsWithCiteseq = renameGenesInSeuratADT(msFNAmerged.samsWithCiteseq,"PTPRC","CD45RA")
msFNAmerged.samsWithCiteseq = renameGenesInSeuratADT(msFNAmerged.samsWithCiteseq,"PTPRC.1","CD45RO")
msFNAmerged.samsWithCiteseq = renameGenesInSeuratADT(msFNAmerged.samsWithCiteseq,"PTPRC.2","CD45")

#rownames(msFNAmerged.samsWithCiteseq@assays[["ADT"]])


# remove all cells from cluster 32, they are erythrocytes and we decided to remove them. 
msFNAmerged <- subset(msFNAmerged,subset = rnaclustersvsAdtManAnn2MaxAnn2 != "Erythroblast")
msFNAmerged.goodSamples <- subset(msFNAmerged.goodSamples,subset = rnaclustersvsAdtManAnn2MaxAnn2 != "Erythroblast")
msFNAmerged.samsWithCiteseq <- subset(msFNAmerged.samsWithCiteseq,subset = rnaclustersvsAdtManAnn2MaxAnn2 != "Erythroblast")


# Renaming for msFNA.merged.tcr and msFNA.merged.bcr has been done in their respective files


#### Save main seurat objects (These rds files can be obtained from EGA) ####
saveRDS(msFNAmerged, file = "data/msFNAmerged.rds")
saveRDS(msFNAmerged.goodSamples, file = "data/msFNAmerged.goodSamples.rds")

saveRDS(msFNAmerged.samsWithCiteseq, file = "data/msFNAmerged.samsWithCiteseq.rds")





















#### Analyses on msFNAmerged.goodSamples ####

DefaultAssay(msFNAmerged) <- "RNA"
DefaultAssay(msFNAmerged.goodSamples) <- "RNA"

## Marker detection for clusters after the final annotation ##

Idents(object = msFNAmerged.goodSamples) <- "rnaclustersvsAdtManAnn2MaxAnn2"
DefaultAssay(msFNAmerged.goodSamples) <- "RNA"

msFNAmerged.FinalAnnotation.RNAdata.markers <- FindAllMarkers(msFNAmerged.goodSamples, only.pos = F, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(msFNAmerged.FinalAnnotation.RNAdata.markers,file="msFNAmerged.FinalAnnotation.RNAdata.markers.csv")

msFNAmerged.FinalAnnotation.RNAdata.markers.top10 <- msFNAmerged.FinalAnnotation.RNAdata.markers %>% group_by(cluster) %>% slice_max(n = 10, order_by = avg_log2FC)
write.csv(msFNAmerged.FinalAnnotation.RNAdata.markers.top10,file="msFNAmerged.FinalAnnotation.RNAdata.markers.top10.csv")


## Pseudobulk DE analysis for each cell type between MS and CTRLs ##

runDEanalysis <- function(seuratObj,countData=NULL,annNames =c("immMain","immFine","rnaClusterAnnotationFromAdt","rnaclustersvsAdtManAnn2MaxAnn2"),excludeBadSamples=F,prefilter=T,resDir = "NULL"){
  
  if(is.null(countData)){
    originalExprData <- GetAssayData(object = seuratObj,assay="RNA",slot = "counts")
  }else{
    originalExprData <- countData
  }
  
  if(!is.null(resDir)){
    dir.create(resDir)
  }
  
  #originalExprData <- rawCountsData
  samples = unique(seuratObj@meta.data$sampleName)
  annotationForDeseq = data.frame(treatment=factor(sapply(samples, function(x) ifelse(grepl("MS",x),"MS","CTRL"))),row.names = samples)
  
  for (cellAnnotType in annNames){
    
    for(unqType in unique(seuratObj[[cellAnnotType]][,1])){
      cellsToAdd = rownames(seuratObj@meta.data[seuratObj[[cellAnnotType]][,1]==unqType,])
      tempMtx <- sapply(samples,function(y) rowSums(originalExprData[,cellsToAdd[grepl(y,cellsToAdd)],drop=F]))
      
      tempMtx <- tempMtx + 1  #to avoid issues in DESeq as a result of some genes with no expression.
      
      if(prefilter==T){
        tempMtx <- tempMtx[rowSums(tempMtx >= 10) >= 3,] # at least X samples with a count of 10 or more, where X can be chosen as the sample size of the smallest group of samples
      }
      
      
      print(unqType)
      
      # exclude bad samples
      if(excludeBadSamples==T){
        tempMtx <- tempMtx[,!colnames(tempMtx) %in% c("MS001","CTRL012","CTRL013","CTRL014")]
        annotationForDeseq <- annotationForDeseq[!rownames(annotationForDeseq) %in% c("MS001","CTRL012","CTRL013","CTRL014"),,drop=F]
      }
      
      dds <- DESeq2::DESeqDataSetFromMatrix(tempMtx, 
                                            colData = annotationForDeseq, 
                                            design = ~ treatment)
      
      dds <- tryCatch(
        expr = {
          DESeq2::DESeq(dds)
        },
        error = function(e){
          message('Using gene-wise dispersion estimate as final estimates!')
          dds <- DESeq2::estimateSizeFactors(dds)
          dds <- DESeq2::estimateDispersionsGeneEst(dds)
          dispersions(dds) <- mcols(dds)$dispGeneEst
          dds <- nbinomWaldTest(dds)
          return(dds)
        }
      )
      
      
      
      contrast <- c("treatment", levels(annotationForDeseq$treatment)[2], levels(annotationForDeseq$treatment)[1])
      
      res <- DESeq2::results(dds, 
                             contrast = contrast,
                             alpha = 0.05,cooksCutoff=FALSE, independentFiltering=FALSE)
      
      
      tempMtxNormFullData = cbind(tempMtx,res)
      tempMtxNormFullData = tempMtxNormFullData[order(tempMtxNormFullData$padj,decreasing=F),]
      
      tempMtxNormFullDataDE <- tempMtxNormFullData[tempMtxNormFullData$padj < 0.05,]
      
      if(grepl("/",unqType)) unqType = gsub("/","-",unqType)
      
      write.csv(tempMtxNormFullData,file=paste0(resDir,"/",cellAnnotType,".",unqType,".pseudobulkDE.csv"))
      if(nrow(tempMtxNormFullDataDE) > 0){
        write.csv(tempMtxNormFullDataDE,file=paste0(resDir,"/",cellAnnotType,".",unqType,".pseudobulkDE_sig.csv"))
      }
    }
    
    
  }
  
}


# DE analysis for our final annotation 
runDEanalysis(seuratObj=msFNAmerged,annNames=c("rnaclustersvsAdtManAnn2MaxAnn2"),resDir = "PseudobulkDE_msFNAmerged_rnaclustersvsAdtManAnn2MaxAnn2")

# DE analysis for our final annotation, excluding bad quality samples
runDEanalysis(seuratObj=msFNAmerged.goodSamples,annNames=c("rnaclustersvsAdtManAnn2MaxAnn2"),resDir = "PseudobulkDE_msFNAmerged.goodSamples_rnaclustersvsAdtManAnn2MaxAnn2")


# plots for the pseudobulk DE results

writeRes="PseudobulkDE_msFNAmerged.goodSamples_rnaclustersvsAdtManAnn2MaxAnn2/"

for(cellT in unique(msFNAmerged.goodSamples@meta.data$rnaclustersvsAdtManAnn2MaxAnn2)){
  
  print(cellT)
  msFNAmerged.goodSamples_ct <- subset(x = msFNAmerged.goodSamples, subset = rnaclustersvsAdtManAnn2MaxAnn2 == cellT)
  
  DEfileName = paste0(writeRes,"rnaclustersvsAdtManAnn2MaxAnn2.",cellT,".pseudobulkDE_sig.csv")
  
  if(file.exists(DEfileName)){
    cellTDE = read.csv(paste0(writeRes,"rnaclustersvsAdtManAnn2MaxAnn2.",cellT,".pseudobulkDE.csv"),row.names = 1)
    cellTDE_sig =  read.csv(paste0(writeRes,"rnaclustersvsAdtManAnn2MaxAnn2.",cellT,".pseudobulkDE_sig.csv"),row.names = 1)
    cellTDE_sig <- cellTDE_sig[order(cellTDE_sig$log2FoldChange,decreasing=T),]
    
    if(nrow(cellTDE_sig) > 0){
      
      plotsDir = paste0(writeRes,cellT,"_plots/")
      dir.create(plotsDir)
      
      EnhancedVolcano(cellTDE,
                      lab = rownames(cellTDE),
                      title = NULL,
                      x = 'log2FoldChange',
                      y = 'padj', pCutoff = 0.05,FCcutoff = 1, maxoverlapsConnectors = Inf, 
                      pointSize = 0.5,
                      labSize = 3,
                      labCol = 'black',
                      labFace = 'bold',
                      boxedLabels = T,
                      colAlpha = 1,
                      legendIconSize = 4.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5) 
      ggsave(paste0(plotsDir,cellT,"_pseudobulkDE_volcanoPlot.pdf"), height = 210, width = 300, units = "mm")
      
      DeGenes = c(rownames(head(cellTDE_sig)),rownames(tail(cellTDE_sig)))
      
      for(degene in DeGenes){
        #degene = "AHNAK"
        
        p1 <- VlnPlot(msFNAmerged.goodSamples_ct, features = degene,group.by = "sampleGroup")
        p2 <- VlnPlot(msFNAmerged.goodSamples_ct, features = degene,group.by = "sampleName")
        p3 <- FeaturePlot(msFNAmerged.goodSamples_ct, features = degene,split.by = "sampleGroup")
        p4 <- RidgePlot(msFNAmerged.goodSamples_ct, features = degene,group.by = "sampleGroup")
        
        cowplot::plot_grid(p1, p2,p3,p4, labels = "AUTO")
        ggsave(paste0(plotsDir,cellT,"_pseudobulkDE_",degene,".pdf"), height = 210, width = 420, units = "mm")
        
      }
      
    }
    
  }
  
  
}



## single cell level DE analysis ##

writeRes="singleCellDE_msFNAmerged.goodSamples_rnaclustersvsAdtManAnn2MaxAnn2/"

for(cellT in unique(msFNAmerged.goodSamples@meta.data$rnaclustersvsAdtManAnn2MaxAnn2)){
  
  print(cellT)
  msFNAmerged.goodSamples_ct <- subset(x = msFNAmerged.goodSamples, subset = rnaclustersvsAdtManAnn2MaxAnn2 == cellT)
  msFNAmerged.samsWithCiteseq_ct <- subset(x = msFNAmerged.samsWithCiteseq, subset = rnaclustersvsAdtManAnn2MaxAnn2 == cellT)
  DefaultAssay(msFNAmerged.samsWithCiteseq_ct) <- "ADT"
  
  # single cell DE
  cellTDE = FindMarkers(msFNAmerged.goodSamples_ct,ident.1="MS",ident.2="CTRL",group.by = "sampleGroup",test.use="MAST")
  cellTDE_sig = cellTDE[cellTDE$p_val_adj < 0.05,]
  cellTDE_sig <- cellTDE_sig[order(cellTDE_sig$avg_log2FC,decreasing=T),]
  
  write.csv(cellTDE,file=paste0(writeRes,"rnaclustersvsAdtManAnn2MaxAnn2.",cellT,".singleCellDE.csv"))
  write.csv(cellTDE_sig,file=paste0(writeRes,"rnaclustersvsAdtManAnn2MaxAnn2.",cellT,".singleCellDE_sig.csv"))
  
  # single cell DE ADT
  ADT_cellTDE = try(FindMarkers(msFNAmerged.samsWithCiteseq_ct,ident.1="MS",ident.2="CTRL",group.by = "sampleGroup",test.use="MAST"),
                    silent = T)
  if(!is.null(nrow(ADT_cellTDE))){
    ADT_cellTDE_sig = ADT_cellTDE[ADT_cellTDE$p_val_adj < 0.05,]
    ADT_cellTDE_sig <- ADT_cellTDE_sig[order(ADT_cellTDE_sig$avg_log2FC,decreasing=T),]
    
    write.csv(ADT_cellTDE,file=paste0(writeRes,"ADT_rnaclustersvsAdtManAnn2MaxAnn2.",cellT,".singleCellDE.csv"))
    write.csv(ADT_cellTDE_sig,file=paste0(writeRes,"ADT_rnaclustersvsAdtManAnn2MaxAnn2.",cellT,".singleCellDE_sig.csv"))
  }else{
    ADT_cellTDE_sig = NULL
  }
  
  # violin plots
  
  if(nrow(cellTDE_sig) > 0){
    
    plotsDir = paste0(writeRes,cellT,"_plots/")
    dir.create(plotsDir)
    
    EnhancedVolcano(cellTDE,
                    lab = rownames(cellTDE),
                    title = NULL,
                    x = 'avg_log2FC',
                    y = 'p_val_adj', pCutoff = 0.05,FCcutoff = 1, maxoverlapsConnectors = Inf, 
                    pointSize = 0.5,
                    labSize = 3,
                    labCol = 'black',
                    labFace = 'bold',
                    boxedLabels = T,
                    colAlpha = 1,
                    legendIconSize = 4.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.5) 
    ggsave(paste0(plotsDir,cellT,"_singleCellDE_volcanoPlot.pdf"), height = 210, width = 300, units = "mm")
    
    DeGenes = c(rownames(head(cellTDE_sig)),rownames(tail(cellTDE_sig)))
    
    for(degene in DeGenes){
      #degene = "LGALS1"
      
      p1 <- VlnPlot(msFNAmerged.goodSamples_ct, features = degene,group.by = "sampleGroup")
      p2 <- VlnPlot(msFNAmerged.goodSamples_ct, features = degene,group.by = "sampleName")
      p3 <- FeaturePlot(msFNAmerged.goodSamples_ct, features = degene,split.by = "sampleGroup")
      p4 <- RidgePlot(msFNAmerged.goodSamples_ct, features = degene,group.by = "sampleGroup")
      
      cowplot::plot_grid(p1, p2,p3,p4, labels = "AUTO")
      ggsave(paste0(plotsDir,cellT,"_singleCellDE_",degene,".pdf"), height = 210, width = 420, units = "mm")
      
    }
    
  }
  
  if(! is.null(ADT_cellTDE_sig)){
    if(nrow(ADT_cellTDE_sig) > 0 ){
      plotsDir = paste0(writeRes,cellT,"_plots/")
      dir.create(plotsDir)
      
      EnhancedVolcano(ADT_cellTDE,
                      lab = rownames(ADT_cellTDE),
                      title = NULL,
                      x = 'avg_log2FC',
                      y = 'p_val_adj', pCutoff = 0.05,FCcutoff = 0.4, maxoverlapsConnectors = Inf, 
                      pointSize = 0.5,
                      labSize = 3,
                      labCol = 'black',
                      labFace = 'bold',
                      boxedLabels = T,
                      colAlpha = 1,
                      legendIconSize = 4.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5) 
      ggsave(paste0(plotsDir,cellT,"_ADT_singleCellDE_volcanoPlot.pdf"), height = 210, width = 300, units = "mm")
      
      DefaultAssay(msFNAmerged.samsWithCiteseq_ct) <- "ADT"
      DoHeatmap(subset(msFNAmerged.samsWithCiteseq_ct, downsample = 100), features = rownames(ADT_cellTDE_sig), group.by="sampleGroup", size = 5)
      ggsave(paste0(plotsDir,cellT,"_ADT_singleCellDE_Heatmap.pdf"), height = 210, width = 300, units = "mm")
      
      
      DeGenesADT = c(rownames(head(ADT_cellTDE_sig)),rownames(tail(ADT_cellTDE_sig)))
      
      for(degeneADT in DeGenesADT){
        #degene = "LGALS1"
        
        p1 <- VlnPlot(msFNAmerged.samsWithCiteseq_ct, assay="ADT",features = degeneADT,group.by = "sampleGroup")
        p2 <- VlnPlot(msFNAmerged.samsWithCiteseq_ct, assay="ADT",features = degeneADT,group.by = "sampleName")
        
        #p3 <- FeaturePlot(msFNAmerged.samsWithCiteseq_ct, features = degeneADT,split.by = "sampleGroup")
        p4 <- RidgePlot(msFNAmerged.samsWithCiteseq_ct,assay="ADT", features = degeneADT,group.by = "sampleGroup")
        
        
        cowplot::plot_grid(p1, p2,p4, labels = "AUTO")
        ggsave(paste0(plotsDir,cellT,"-ADT_singleCellDE_",degeneADT,".pdf"), height = 210, width = 420, units = "mm")
        
      }
    }
  }
  
  
}  




## Detect single cell DE genes between MS vs CTRLs for for B cell types ## 

# For GC B cells sub types
# Seurat object msFNAmerged.goodSamples.bcells4 contains subclustering of GC B cells after they were subset from the full data

msFNAmerged.goodSamples.bcells4.js = readRDS("data/msFNAmerged.goodSamples.bcells4")

# For Subclusters of GC B cells
writeRes="singleCellDE_GCBcell_subclusters/"
dir.create(writeRes)

for(cellT in unique(msFNAmerged.goodSamples.bcells4.js@meta.data$GC_type)){
  
  print(cellT)
  msFNAmerged.goodSamples.bcells4_ct <- subset(x = msFNAmerged.goodSamples.bcells4.js, subset = GC_type == cellT)
  msFNAmerged.samsWithCiteseq.bcells4_ct <- msFNAmerged.samsWithCiteseq[,rownames(msFNAmerged.goodSamples.bcells4_ct@meta.data)]
  DefaultAssay(msFNAmerged.samsWithCiteseq.bcells4_ct) <- "ADT"
  
  
  
  # single cell DE
  cellTDE = FindMarkers(msFNAmerged.goodSamples.bcells4_ct,ident.1="MS",ident.2="CTRL",group.by = "sampleGroup",test.use="MAST")
  cellTDE_sig = cellTDE[cellTDE$p_val_adj < 0.05,]
  cellTDE_sig <- cellTDE_sig[order(cellTDE_sig$avg_log2FC,decreasing=T),]
  
  write.csv(cellTDE,file=paste0(writeRes,"GCbcell.subcluster.",cellT,".singleCellDE.csv"))
  write.csv(cellTDE_sig,file=paste0(writeRes,"GCbcell.subcluster.",cellT,".singleCellDE_sig.csv"))
  
  # single cell DE ADT
  ADT_cellTDE = try(FindMarkers(msFNAmerged.samsWithCiteseq.bcells4_ct,ident.1="MS",ident.2="CTRL",group.by = "sampleGroup",test.use="MAST"),
                    silent = T)
  if(!is.null(nrow(ADT_cellTDE))){
    ADT_cellTDE_sig = ADT_cellTDE[ADT_cellTDE$p_val_adj < 0.05,]
    ADT_cellTDE_sig <- ADT_cellTDE_sig[order(ADT_cellTDE_sig$avg_log2FC,decreasing=T),]
    
    write.csv(ADT_cellTDE,file=paste0(writeRes,"ADT_GCbcell.subcluster.",cellT,".singleCellDE.csv"))
    write.csv(ADT_cellTDE_sig,file=paste0(writeRes,"ADT_GCbcell.subcluster.",cellT,".singleCellDE_sig.csv"))
  }else{
    ADT_cellTDE_sig = NULL
  }
  
  # violin plots
  
  if(nrow(cellTDE_sig) > 0){
    
    plotsDir = paste0(writeRes,cellT,"_plots/")
    dir.create(plotsDir)
    
    EnhancedVolcano(cellTDE,
                    lab = rownames(cellTDE),
                    title = NULL,
                    x = 'avg_log2FC',
                    y = 'p_val_adj', pCutoff = 0.05,FCcutoff = 1, maxoverlapsConnectors = Inf, 
                    pointSize = 0.5,
                    labSize = 3,
                    labCol = 'black',
                    labFace = 'bold',
                    boxedLabels = T,
                    colAlpha = 1,
                    legendIconSize = 4.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.5) 
    ggsave(paste0(plotsDir,cellT,"_singleCellDE_volcanoPlot.pdf"), height = 210, width = 300, units = "mm")
    
    DeGenes = c(rownames(head(cellTDE_sig)),rownames(tail(cellTDE_sig)))
    
    for(degene in DeGenes){
      #degene = "LGALS1"
      
      p1 <- VlnPlot(msFNAmerged.goodSamples.bcells4_ct, features = degene,group.by = "sampleGroup")
      p2 <- VlnPlot(msFNAmerged.goodSamples.bcells4_ct, features = degene,group.by = "sampleName")
      p3 <- FeaturePlot(msFNAmerged.goodSamples.bcells4_ct, features = degene,split.by = "sampleGroup")
      p4 <- RidgePlot(msFNAmerged.goodSamples.bcells4_ct, features = degene,group.by = "sampleGroup")
      
      cowplot::plot_grid(p1, p2,p3,p4, labels = "AUTO")
      ggsave(paste0(plotsDir,cellT,"_singleCellDE_",degene,".pdf"), height = 210, width = 420, units = "mm")
      
    }
    
  }
  
  if(! is.null(ADT_cellTDE_sig)){
    if(nrow(ADT_cellTDE_sig) > 0 ){
      plotsDir = paste0(writeRes,cellT,"_plots/")
      dir.create(plotsDir)
      
      EnhancedVolcano(ADT_cellTDE,
                      lab = rownames(ADT_cellTDE),
                      title = NULL,
                      x = 'avg_log2FC',
                      y = 'p_val_adj', pCutoff = 0.05,FCcutoff = 0.4, maxoverlapsConnectors = Inf, 
                      pointSize = 0.5,
                      labSize = 3,
                      labCol = 'black',
                      labFace = 'bold',
                      boxedLabels = T,
                      colAlpha = 1,
                      legendIconSize = 4.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5) 
      ggsave(paste0(plotsDir,cellT,"_ADT_singleCellDE_volcanoPlot.pdf"), height = 210, width = 300, units = "mm")
      
      DefaultAssay(msFNAmerged.samsWithCiteseq.bcells4_ct) <- "ADT"
      DoHeatmap(subset(msFNAmerged.samsWithCiteseq.bcells4_ct, downsample = 100), features = rownames(ADT_cellTDE_sig), group.by="sampleGroup", size = 5)
      ggsave(paste0(plotsDir,cellT,"_ADT_singleCellDE_Heatmap.pdf"), height = 210, width = 300, units = "mm")
      
      
      DeGenesADT = c(rownames(head(ADT_cellTDE_sig)),rownames(tail(ADT_cellTDE_sig)))
      
      for(degeneADT in DeGenesADT){
        #degene = "LGALS1"
        
        p1 <- VlnPlot(msFNAmerged.samsWithCiteseq.bcells4_ct, assay="ADT",features = degeneADT,group.by = "sampleGroup")
        p2 <- VlnPlot(msFNAmerged.samsWithCiteseq.bcells4_ct, assay="ADT",features = degeneADT,group.by = "sampleName")
        
        #p3 <- FeaturePlot(msFNAmerged.samsWithCiteseq_ct, features = degeneADT,split.by = "sampleGroup")
        p4 <- RidgePlot(msFNAmerged.samsWithCiteseq.bcells4_ct,assay="ADT", features = degeneADT,group.by = "sampleGroup")
        
        
        cowplot::plot_grid(p1, p2,p4, labels = "AUTO")
        ggsave(paste0(plotsDir,cellT,"-ADT_singleCellDE_",degeneADT,".pdf"), height = 210, width = 420, units = "mm")
        
      }
    }
  }
  
  
}  


# For mem and naive B cells
# Seurat object msFNAmerged.goodSamples.membcells contains subclustering of naive and memory B cells after they were subset from the full data

msFNAmerged.goodSamples.memNaivebcells.js = readRDS("data/msFNAmerged.goodSamples.membcells")

writeRes="singleCellDE_memorynaiveB_subclusters/"
dir.create(writeRes)

for(cellT in unique(msFNAmerged.goodSamples.memNaivebcells.js@meta.data$memnaive_type)){
  
  print(cellT)
  msFNAmerged.goodSamples.bcells4_ct <- subset(x = msFNAmerged.goodSamples.memNaivebcells.js, subset = memnaive_type == cellT)
  msFNAmerged.samsWithCiteseq.bcells4_ct <- msFNAmerged.samsWithCiteseq[,rownames(msFNAmerged.goodSamples.bcells4_ct@meta.data)]
  DefaultAssay(msFNAmerged.samsWithCiteseq.bcells4_ct) <- "ADT"
  
  
  
  # single cell DE
  cellTDE = FindMarkers(msFNAmerged.goodSamples.bcells4_ct,ident.1="MS",ident.2="CTRL",group.by = "sampleGroup",test.use="MAST")
  cellTDE_sig = cellTDE[cellTDE$p_val_adj < 0.05,]
  cellTDE_sig <- cellTDE_sig[order(cellTDE_sig$avg_log2FC,decreasing=T),]
  
  write.csv(cellTDE,file=paste0(writeRes,"GCbcell.subcluster.",cellT,".singleCellDE.csv"))
  write.csv(cellTDE_sig,file=paste0(writeRes,"GCbcell.subcluster.",cellT,".singleCellDE_sig.csv"))
  
  # single cell DE ADT
  ADT_cellTDE = try(FindMarkers(msFNAmerged.samsWithCiteseq.bcells4_ct,ident.1="MS",ident.2="CTRL",group.by = "sampleGroup",test.use="MAST"),
                    silent = T)
  if(!is.null(nrow(ADT_cellTDE))){
    ADT_cellTDE_sig = ADT_cellTDE[ADT_cellTDE$p_val_adj < 0.05,]
    ADT_cellTDE_sig <- ADT_cellTDE_sig[order(ADT_cellTDE_sig$avg_log2FC,decreasing=T),]
    
    write.csv(ADT_cellTDE,file=paste0(writeRes,"ADT_GCbcell.subcluster.",cellT,".singleCellDE.csv"))
    write.csv(ADT_cellTDE_sig,file=paste0(writeRes,"ADT_GCbcell.subcluster.",cellT,".singleCellDE_sig.csv"))
  }else{
    ADT_cellTDE_sig = NULL
  }
  
  # violin plots
  
  if(nrow(cellTDE_sig) > 0){
    
    plotsDir = paste0(writeRes,cellT,"_plots/")
    dir.create(plotsDir)
    
    EnhancedVolcano(cellTDE,
                    lab = rownames(cellTDE),
                    title = NULL,
                    x = 'avg_log2FC',
                    y = 'p_val_adj', pCutoff = 0.05,FCcutoff = 1, maxoverlapsConnectors = Inf, 
                    pointSize = 0.5,
                    labSize = 3,
                    labCol = 'black',
                    labFace = 'bold',
                    boxedLabels = T,
                    colAlpha = 1,
                    legendIconSize = 4.0,
                    drawConnectors = TRUE,
                    widthConnectors = 0.5) 
    ggsave(paste0(plotsDir,cellT,"_singleCellDE_volcanoPlot.pdf"), height = 210, width = 300, units = "mm")
    
    DeGenes = c(rownames(head(cellTDE_sig)),rownames(tail(cellTDE_sig)))
    
    for(degene in DeGenes){
      #degene = "LGALS1"
      
      p1 <- VlnPlot(msFNAmerged.goodSamples.bcells4_ct, features = degene,group.by = "sampleGroup")
      p2 <- VlnPlot(msFNAmerged.goodSamples.bcells4_ct, features = degene,group.by = "sampleName")
      p3 <- FeaturePlot(msFNAmerged.goodSamples.bcells4_ct, features = degene,split.by = "sampleGroup")
      p4 <- RidgePlot(msFNAmerged.goodSamples.bcells4_ct, features = degene,group.by = "sampleGroup")
      
      cowplot::plot_grid(p1, p2,p3,p4, labels = "AUTO")
      ggsave(paste0(plotsDir,cellT,"_singleCellDE_",degene,".pdf"), height = 210, width = 420, units = "mm")
      
    }
    
  }
  
  if(! is.null(ADT_cellTDE_sig)){
    if(nrow(ADT_cellTDE_sig) > 0 ){
      plotsDir = paste0(writeRes,cellT,"_plots/")
      dir.create(plotsDir)
      
      EnhancedVolcano(ADT_cellTDE,
                      lab = rownames(ADT_cellTDE),
                      title = NULL,
                      x = 'avg_log2FC',
                      y = 'p_val_adj', pCutoff = 0.05,FCcutoff = 0.4, maxoverlapsConnectors = Inf, 
                      pointSize = 0.5,
                      labSize = 3,
                      labCol = 'black',
                      labFace = 'bold',
                      boxedLabels = T,
                      colAlpha = 1,
                      legendIconSize = 4.0,
                      drawConnectors = TRUE,
                      widthConnectors = 0.5) 
      ggsave(paste0(plotsDir,cellT,"_ADT_singleCellDE_volcanoPlot.pdf"), height = 210, width = 300, units = "mm")
      
      DefaultAssay(msFNAmerged.samsWithCiteseq.bcells4_ct) <- "ADT"
      DoHeatmap(subset(msFNAmerged.samsWithCiteseq.bcells4_ct, downsample = 100), features = rownames(ADT_cellTDE_sig), group.by="sampleGroup", size = 5)
      ggsave(paste0(plotsDir,cellT,"_ADT_singleCellDE_Heatmap.pdf"), height = 210, width = 300, units = "mm")
      
      
      DeGenesADT = c(rownames(head(ADT_cellTDE_sig)),rownames(tail(ADT_cellTDE_sig)))
      
      for(degeneADT in DeGenesADT){
        #degene = "LGALS1"
        
        p1 <- VlnPlot(msFNAmerged.samsWithCiteseq.bcells4_ct, assay="ADT",features = degeneADT,group.by = "sampleGroup")
        p2 <- VlnPlot(msFNAmerged.samsWithCiteseq.bcells4_ct, assay="ADT",features = degeneADT,group.by = "sampleName")
        
        #p3 <- FeaturePlot(msFNAmerged.samsWithCiteseq_ct, features = degeneADT,split.by = "sampleGroup")
        p4 <- RidgePlot(msFNAmerged.samsWithCiteseq.bcells4_ct,assay="ADT", features = degeneADT,group.by = "sampleGroup")
        
        
        cowplot::plot_grid(p1, p2,p4, labels = "AUTO")
        ggsave(paste0(plotsDir,cellT,"-ADT_singleCellDE_",degeneADT,".pdf"), height = 210, width = 420, units = "mm")
        
      }
    }
  }
  
  
}  




## Compare differential cell type abundance ##

annNames = c("rnaclustersvsAdtManAnn2MaxAnn2","immFine")

for(cellAnnotType in annNames){
  
  cellTypeVersusSampleCount = table(msFNAmerged.goodSamples[[cellAnnotType]][,1],msFNAmerged.goodSamples@meta.data$sampleName)
  ctrlSamples <- grepl("CTRL",colnames(cellTypeVersusSampleCount))
  msSamples <- !grepl("CTRL",colnames(cellTypeVersusSampleCount))
  
  cellTypeVersusSampleCountNorm <- apply(cellTypeVersusSampleCount + 1,2, function(x) (x/sum(x)))
  
  
  # perform unpaired t.test
  diffCellAbundance = t(apply(cellTypeVersusSampleCountNorm,1,function(x) c(t.test(x[msSamples],x[ctrlSamples])$p.value,
                                                                            log2(mean(x[msSamples])/mean(x[ctrlSamples])),lsr::cohensD(x[msSamples]/x[ctrlSamples])
                                                                            )))
  
  
  colnames(diffCellAbundance) <- c("unpaired two-samples t.test.p.value","log2FC(meanMS/meanCTRL)","EffectSize(Cohens_d)")
  
  fullDiffCellData = cbind(cellTypeVersusSampleCount,cellTypeVersusSampleCountNorm,diffCellAbundance)
  fullDiffCellData <- fullDiffCellData[order(fullDiffCellData[,which(grepl("p.value",colnames(fullDiffCellData)))],decreasing=F),]
  
  
  write.csv(fullDiffCellData,paste0("mergedDataDiferentialCellTypeAbundance_",cellAnnotType,"_ttest.differentialCellTypeAbudance.csv"))
  
  
}



## Plots ##

# Plot 1: Landscape of immune cells in the cervical lymph nodes

DimPlot(msFNAmerged.goodSamples,label=T,group.by = "sampleGroup",repel =T)
ggsave("UmapPlot_msFNAmerged.goodSamples.sampleGroups.pdf",height = 210, width = 320, units = "mm")


DimPlot(msFNAmerged.goodSamples,label=T,group.by = "rnaclustersvsAdtManAnn2MaxAnn2",label.size = 3, split.by = "sampleGroup",
        cols=DiscretePalette(19, palette = "polychrome"))
ggsave("UmapPlot_msFNAmerged.goodSamples.AnnotationsBySampleGroups.pdf",height = 210, width = 320, units = "mm")



DimPlot(msFNAmerged.goodSamples,label=T,group.by = "rnaclustersvsAdtManAnn2MaxAnn2",label.size = 4,
        cols=DiscretePalette(19, palette = "polychrome"))
ggsave("UmapPlot_msFNAmerged.goodSamples.Annotations.pdf",height = 210, width = 320, units = "mm")



# Plot2: Bar plot of the overall frequency of cell types

pdf(file="msFNAmerged.goodSamples_CellTypesCounts.pdf")
par(mar = c(9, 4.1, 4.1, 2.1))
barplot(sort(table(msFNAmerged.goodSamples@meta.data$rnaclustersvsAdtManAnn2MaxAnn2),decreasing=T),
        las = 2,
        ylim=c(0,30000),ylab="Detected cell types")

dev.off()

#--
pdf(file="msFNAmerged.goodSamples_CellTypeProportions.pdf")
par(mar = c(9, 4.1, 4.1, 2.1))
barplot(100 * (sort(table(msFNAmerged.goodSamples@meta.data$rnaclustersvsAdtManAnn2MaxAnn2),decreasing=T)/sum(table(msFNAmerged.goodSamples@meta.data$rnaclustersvsAdtManAnn2MaxAnn2))),
        las = 2,
        ylim=c(0,40),ylab="Detected cell types (%)")

dev.off()


#-- controls and MS separately
pdf(file="msFNAmerged.goodSamples_CellTypeProportionsBySampleGroupsMS.pdf")
par(mar = c(9, 4.1, 4.1, 2.1))


barplot(100 * (sort(table(msFNAmerged.goodSamples@meta.data[msFNAmerged.goodSamples@meta.data$sampleGroup=="MS",]$rnaclustersvsAdtManAnn2MaxAnn2),decreasing=T)/sum(table(msFNAmerged.goodSamples@meta.data[msFNAmerged.goodSamples@meta.data$sampleGroup=="MS",]$rnaclustersvsAdtManAnn2MaxAnn2))),
        las = 2,
        ylim=c(0,40),ylab="Detected cell types (%)")

dev.off()


pdf(file="msFNAmerged.goodSamples_CellTypeProportionsBySampleGroupsCTRL.pdf")
par(mar = c(9, 4.1, 4.1, 2.1))
barplot(100 * (sort(table(msFNAmerged.goodSamples@meta.data[msFNAmerged.goodSamples@meta.data$sampleGroup=="CTRL",]$rnaclustersvsAdtManAnn2MaxAnn2),decreasing=T)/sum(table(msFNAmerged.goodSamples@meta.data[msFNAmerged.goodSamples@meta.data$sampleGroup=="CTRL",]$rnaclustersvsAdtManAnn2MaxAnn2))),
        las = 2,
        ylim=c(0,50),ylab="Detected cell types (%)")

dev.off()



##1. Compare frequency distribution of the cell types between the two groups (MS011 included).
cellCountsInMS = table(msFNAmerged.goodSamples@meta.data[msFNAmerged.goodSamples@meta.data$sampleGroup=="MS",]$rnaclustersvsAdtManAnn2MaxAnn2)
cellCountsInCTRL = table(msFNAmerged.goodSamples@meta.data[msFNAmerged.goodSamples@meta.data$sampleGroup=="CTRL",]$rnaclustersvsAdtManAnn2MaxAnn2)

cellFrequenyData = as.matrix(cbind(MS=cellCountsInMS,CTRL=cellCountsInCTRL))
test <- chisq.test(cellFrequenyData)

M2 <- as.table(rbind(cellCountsInMS, cellCountsInCTRL))
dimnames(M2) <- list(sampleGroup = c("MS", "CTRL"),
                     CellType = names(cellCountsInMS))
(Xsq2 <- chisq.test(t(M2))) 
Xsq2$observed  
Xsq2$expected   
Xsq2$residuals  
Xsq2$stdres 

stdResidualSig_cutoff <- function(nrow = 2, ncol = 2, alpha = .01){
  value <- qnorm(p=1-((alpha/2)/(nrow*ncol)))
  #output <- paste("Residual cutoff is", round(value, 3))
  return(round(value, 3))
}


pdf(file="msFNAmerged.goodSamples_cellTypePropCombinedDF_chisquare_standardizedResiduals.pdf")
par(mar = c(9, 4.1, 4.1, 2.1))

barplot(sort(Xsq2$residuals[,"MS"]), las = 2,  ylim=c(-15,15),
        ylab="Standardized residual (MS)")
abline(h=c(stdResidualSig_cutoff(2,19),-stdResidualSig_cutoff(2,19)),col="red",lty=2)
axis(side=2, at=c(stdResidualSig_cutoff(2,19),-stdResidualSig_cutoff(2,19)), 
     labels=c(round(stdResidualSig_cutoff(2,19),2),round(-stdResidualSig_cutoff(2,19),2)),las=2,cex.axis=0.5)
dev.off()




# combined plot of cell type proportions per sample group
cellPropInMS = 100 * (sort(table(msFNAmerged.goodSamples@meta.data[msFNAmerged.goodSamples@meta.data$sampleGroup=="MS",]$rnaclustersvsAdtManAnn2MaxAnn2),decreasing=T)/sum(table(msFNAmerged.goodSamples@meta.data[msFNAmerged.goodSamples@meta.data$sampleGroup=="MS",]$rnaclustersvsAdtManAnn2MaxAnn2)))
cellPropInCTRL = 100 * (sort(table(msFNAmerged.goodSamples@meta.data[msFNAmerged.goodSamples@meta.data$sampleGroup=="CTRL",]$rnaclustersvsAdtManAnn2MaxAnn2),decreasing=T)/sum(table(msFNAmerged.goodSamples@meta.data[msFNAmerged.goodSamples@meta.data$sampleGroup=="CTRL",]$rnaclustersvsAdtManAnn2MaxAnn2)))

cellPropInMS_DF = data.frame(cellPropInMS)
colnames(cellPropInMS_DF) <- c("CellType","Relative_Frequency")
cellPropInMS_DF$sampleGroup="MS"

cellPropInCTRL_DF = data.frame(cellPropInCTRL)
colnames(cellPropInCTRL_DF) <- c("CellType","Relative_Frequency")
cellPropInCTRL_DF$sampleGroup="CTRL"

cellTypePropCombinedDF = rbind(cellPropInMS_DF,cellPropInCTRL_DF)



ggplot(data=cellTypePropCombinedDF, aes(x=CellType, y=Relative_Frequency, fill=sampleGroup)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_manual(values=c("lightblue","Brown")) + 
  theme(legend.position = "top") +  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ylab("Relative Frequency of detected cell types (%)") + xlab("")

ggsave("msFNAmerged.goodSamples_cellTypePropCombinedDF.pdf", height = 210, width = 320, units = "mm")



##2. Compare frequency distribution of the cell types between the two groups (MS011 excluded).

## MS011 is in active relapse, and shows distinct characteristics compared to other patients. We wanted to check if the general
## cell type frequency distribution observed in patients remains the same when we remove this patient from the analysis.

msFNAmerged.goodSamples.minusMS011 <- subset(msFNAmerged.goodSamples, subset = sampleName != "MS011")

cellCountsInMS_minusMS011 = table(msFNAmerged.goodSamples.minusMS011@meta.data[msFNAmerged.goodSamples.minusMS011@meta.data$sampleGroup=="MS",]$rnaclustersvsAdtManAnn2MaxAnn2)
cellCountsInCTRL_minusMS011 = table(msFNAmerged.goodSamples.minusMS011@meta.data[msFNAmerged.goodSamples.minusMS011@meta.data$sampleGroup=="CTRL",]$rnaclustersvsAdtManAnn2MaxAnn2)


M2_minusMS011 <- as.table(rbind(cellCountsInMS_minusMS011, cellCountsInCTRL_minusMS011))
dimnames(M2_minusMS011) <- list(sampleGroup = c("MS", "CTRL"),
                                CellType = names(cellCountsInMS))
(Xsq2 <- chisq.test(t(M2_minusMS011)))  
Xsq2$observed   
Xsq2$expected  
Xsq2$residuals  
Xsq2$stdres 

pdf(file="msFNAmerged.goodSamples_minusMS011_cellTypePropCombinedDF_chisquare_standardizedResiduals.pdf")
par(mar = c(9, 4.1, 4.1, 2.1))

barplot(sort(Xsq2$residuals[,"MS"]), las = 2,  ylim=c(-15,15),
        ylab="Standardized residual (MS)")
abline(h=c(stdResidualSig_cutoff(2,19),-stdResidualSig_cutoff(2,19)),col="red",lty=2)
axis(side=2, at=c(stdResidualSig_cutoff(2,19),-stdResidualSig_cutoff(2,19)), 
     labels=c(round(stdResidualSig_cutoff(2,19),2),round(-stdResidualSig_cutoff(2,19),2)),las=2,cex.axis=0.5)
dev.off()


# combined plot of cell type proportions per sample group
cellPropInMS_minusMS011 = 100 * (sort(table(msFNAmerged.goodSamples.minusMS011@meta.data[msFNAmerged.goodSamples.minusMS011@meta.data$sampleGroup=="MS",]$rnaclustersvsAdtManAnn2MaxAnn2),decreasing=T)/sum(table(msFNAmerged.goodSamples.minusMS011@meta.data[msFNAmerged.goodSamples.minusMS011@meta.data$sampleGroup=="MS",]$rnaclustersvsAdtManAnn2MaxAnn2)))
cellPropInCTRL_minusMS011 = 100 * (sort(table(msFNAmerged.goodSamples.minusMS011@meta.data[msFNAmerged.goodSamples.minusMS011@meta.data$sampleGroup=="CTRL",]$rnaclustersvsAdtManAnn2MaxAnn2),decreasing=T)/sum(table(msFNAmerged.goodSamples.minusMS011@meta.data[msFNAmerged.goodSamples.minusMS011@meta.data$sampleGroup=="CTRL",]$rnaclustersvsAdtManAnn2MaxAnn2)))

cellPropInMS_minusMS011_DF = data.frame(cellPropInMS_minusMS011)
colnames(cellPropInMS_minusMS011_DF) <- c("CellType","Relative_Frequency")
cellPropInMS_minusMS011_DF$sampleGroup="MS"

cellPropInCTRL_minusMS011_DF = data.frame(cellPropInCTRL_minusMS011)
colnames(cellPropInCTRL_minusMS011_DF) <- c("CellType","Relative_Frequency")
cellPropInCTRL_minusMS011_DF$sampleGroup="CTRL"

cellTypePropCombinedDF_minusMS011 = rbind(cellPropInMS_minusMS011_DF,cellPropInCTRL_minusMS011_DF)



ggplot(data=cellTypePropCombinedDF_minusMS011, aes(x=CellType, y=Relative_Frequency, fill=sampleGroup)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  scale_fill_manual(values=c("lightblue","Brown")) + 
  theme(legend.position = "top") +  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ylab("Relative Frequency of detected cell types (%)") + xlab("")

ggsave("msFNAmerged.goodSamples_minusMS011_cellTypePropCombinedDF.pdf", height = 210, width = 320, units = "mm")








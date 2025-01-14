#!/usr/bin/env Rscript

# This runs seurat on a single sample whose 10x filtered data matrix directory (after cellranger count or cellranger multi) passed via arg 

library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)

# read in data name
args <- commandArgs(trailingOnly = TRUE)

dataSample = args[1]
samPrefix = basename(dataSample)

# data sample would be something like this:
#dataSample = "/scratch/project_2005392/Joona/MS_FNA/scRNAseq/MS003/"


samp <- Read10X(dataSample)

if(is.list(samp)){
  samp = samp$`Gene Expression`
}
numGenesDetectedPerCell <- apply(samp,2, function(x) sum(x>0))
lowerFeaturesInitialcutoff = min(numGenesDetectedPerCell)

samSeuratObject <- CreateSeuratObject(samp, project = samPrefix, min.cells = 3, min.features = lowerFeaturesInitialcutoff)



# QC: 
samSeuratObject[["percent.mt"]] <- PercentageFeatureSet(samSeuratObject, pattern = "^MT-")
VlnPlot(samSeuratObject, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
ggsave(paste0(samPrefix,"_vln_qc1.pdf"), width = 10, height = 6)


# scatter plots
plot1 <- FeatureScatter(samSeuratObject, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(samSeuratObject, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
plot1 + plot2
ggsave(paste0(samPrefix,"_ftscat_qc2.pdf"), width = 10, height = 6)



#filtering: we find better cutoff for the maximum nFeature_RNA, and percent.mt
  source("/scratch/project_2005392/Pia/up to date R/seuratPublicDataForComparison.R")

  # the pbmc example data in seurat used lower 98% for mt. We use a cut off that retains similar percentile of the data for mt. And 99 percentile lower #genes detected
  #mt_maxQuantile & nfeatures_max obtained from seuratPublicDataForComparison.R
  
  mtPossibleCutoff = as.numeric(round(quantile(samSeuratObject@meta.data$percent.mt,probs = c(mt_maxQuantile))))
  
  maxFeatures= as.numeric(round(quantile(samSeuratObject@meta.data$nFeature_RNA,probs = c(nfeatures_max))))
  
  
samSeuratObject <- subset(samSeuratObject, subset = nFeature_RNA > lowerFeaturesInitialcutoff & nFeature_RNA < maxFeatures & percent.mt < mtPossibleCutoff)


#normalizing data
samSeuratObject <- NormalizeData(samSeuratObject, normalization.method = "LogNormalize", scale.factor = 10000)

# find variable genes
samSeuratObject <- FindVariableFeatures(samSeuratObject, selection.method = "vst", nfeatures = 2000)


#top 10 highly variable genes
top10.varGenes <- head(VariableFeatures(samSeuratObject), 10)
plot1 <- VariableFeaturePlot(samSeuratObject)
LabelPoints(plot = plot1, points = top10.varGenes, repel = TRUE, xnudge = 0, ynudge = 0)
ggsave(paste0(samPrefix,"_hvg.pdf"), width = 10, height = 6)


# Cluster Analysis:

samSeuratObject <- ScaleData(samSeuratObject, features = rownames(samSeuratObject))
samSeuratObject <- RunPCA(samSeuratObject, features = VariableFeatures(object = samSeuratObject))
# Examine and visualize PCA results a few different ways
ElbowPlot(samSeuratObject)
ggsave(paste0(samPrefix,"_elbow.pdf"), height = 10, width = 6)


##cell cycle scoring
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes

samSeuratObject  <- CellCycleScoring(samSeuratObject, s.features = s.genes, g2m.features = g2m.genes, set.ident = F)


# Running a PCA on cell cycle genes, do cells segregate by phase
samSeuratObject <- RunPCA(samSeuratObject, features = c(s.genes, g2m.genes))
DimPlot(samSeuratObject,group.by="Phase")
ggsave(paste0(samPrefix,"_PCAbyCellCycleGenes.pdf"), height = 10, width = 6)



# Determine number of PCs with Jackstraw
samSeuratObject <- RunPCA(samSeuratObject, features = VariableFeatures(object = samSeuratObject))

samSeuratObject <- JackStraw(samSeuratObject, prop.freq=0.1, dims = 50, num.replicate = 100)
samSeuratObject <- ScoreJackStraw(samSeuratObject, dims = 1:50)
JackStrawPlot(samSeuratObject, dims = 1:50)
ggsave(paste0(samPrefix,"_JackstrawPvalues.pdf"), height = 297, width = 210, units = "mm")

score.df <- JS(object = samSeuratObject[['pca']], slot = 'overall')
selectedPCs <- which(score.df[,2] < 0.05)


# Clustering:
samSeuratObject <- FindNeighbors(samSeuratObject, dims = selectedPCs)
samSeuratObject <- FindClusters(samSeuratObject, resolution = 0.6) 

samSeuratObject <- RunUMAP(samSeuratObject, dims = selectedPCs)


DimPlot(samSeuratObject, reduction = "umap", label = T) + ggtitle(samPrefix)
ggsave(paste0(samPrefix,"_UMAP.pdf"), height = 10, width = 10)



#cluster biomarkers
cluster.markers <- c()
cluster.markers.top10 <- c()


cluster.markers <- FindAllMarkers(samSeuratObject, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cluster.markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC) -> cluster.markers.top10

write.csv(cluster.markers,file=paste0(samPrefix,"_clusterMarkers.csv"))
write.csv(cluster.markers.top10,file=paste0(samPrefix,"_top10_clusterMarkers.csv"))


DefaultAssay(samSeuratObject) <- "RNA"




# Automatic cell annotation using singleR and the Monaco reference 

library(SingleR)
library(celldex)
hpca.se <- HumanPrimaryCellAtlasData()


immune.monaco <- MonacoImmuneData()


# annotate cells 

normData <- GetAssayData(samSeuratObject,assay="RNA",slot="data")

#1
pred.hesc <- SingleR(test = normData, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.fine)

samSeuratObject[["SingleR.HumanPrimaryCellAtlasData.fineLabels"]] <- pred.hesc$pruned.labels

#2
pred.hesc <- SingleR(test = normData, ref = hpca.se, assay.type.test=1,
                     labels = hpca.se$label.main)

samSeuratObject[["SingleR.HumanPrimaryCellAtlasData.mainLabels"]] <- pred.hesc$pruned.labels


#3
pred.hesc <- SingleR(test = normData, ref = immune.monaco, assay.type.test=1,
                     labels = immune.monaco$label.fine)

samSeuratObject[["SingleR.immune.monaco.fineLabels"]] <- pred.hesc$pruned.labels

#4
pred.hesc <- SingleR(test = normData, ref = immune.monaco, assay.type.test=1,
                     labels = immune.monaco$label.main)

samSeuratObject[["SingleR.immune.monaco.mainLabels"]] <- pred.hesc$pruned.labels


DimPlot(samSeuratObject, reduction = "umap", label = TRUE, group.by=c("seurat_clusters","SingleR.immune.monaco.mainLabels"), repel=TRUE, pt.size= 0.5,label.size = 5) + NoLegend() + ggtitle(samPrefix)

ggsave(paste0(samPrefix,"_clusterAnnotation.pdf"), height = 210, width = 350, units = "mm")


# save analysis
write.csv(samSeuratObject@meta.data,file=paste0(samPrefix,"_complete_metaData.csv"))


# cell composition of main cells in the sample
cellComposition = as.data.frame(table(samSeuratObject@meta.data$SingleR.HumanPrimaryCellAtlasData.mainLabels))

cellComposition$proportionsInPerc = (cellComposition$Freq/sum(cellComposition$Freq))

write.csv(cellComposition,
          file=paste0(samPrefix,"_cellCompositionHPCA.csv"))


# cell composition of finer annotation using the immune monaco reference
cellCompositionFine = as.data.frame(table(samSeuratObject@meta.data$SingleR.immune.monaco.fineLabels))

cellCompositionFine$proportionsInPerc = (cellCompositionFine$Freq/sum(cellCompositionFine$Freq))

write.csv(cellCompositionFine,
          file=paste0(samPrefix,"_cellCompositionFineImmuneMonaco.csv"))


# save analysis
saveRDS(samSeuratObject, file = paste0(samPrefix,"_seuratAnalysis.rds"))





#packages
library(dplyr)
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(ggpubr)

library(scRepertoire)
source("TCR_functions.R")


##first read in all data
msFNAmerged.goodSamples <- readRDS(file = "data/msFNAmerged.goodSamples.rds")
msFNA.merged.tcr <- readRDS(file = "data/msFNA.merged.tcr.rds")

msFNAmerged.goodSamples <- SetIdent(msFNAmerged.goodSamples, value = "rnaclustersvsAdtManAnn2MaxAnn2")

#### Subcluster the CD8s ####

unique(msFNAmerged.goodSamples$rnaclustersvsAdtManAnn2MaxAnn2)
head(msFNAmerged.goodSamples)
cd8_df <- subset(x = msFNAmerged.goodSamples, subset = rnaclustersvsAdtManAnn2MaxAnn2 %in% c("Naive CD8", "Memory CD8")) #this includes CD8 T cells

cd8_df <- NormalizeData(cd8_df)
cd8_df <- FindVariableFeatures(cd8_df, selection.method = "vst", nfeatures = 3000)
all.genes.cd8_df <- rownames(cd8_df)
cd8_df <- ScaleData(cd8_df, features = all.genes.cd8_df)
top100 <- head(VariableFeatures(cd8_df), 100)
vargenes_cd8 <- VariableFeatures(object = cd8_df)

#cell cycle
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cd8_df <- CellCycleScoring(cd8_df, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

#dimensionality and reduction
cd8_df <- RunPCA(cd8_df, features = VariableFeatures(object = cd8_df))
VizDimLoadings(cd8_df, dims = 1:2, reduction = "pca")
DimPlot(cd8_df, reduction = "pca")
DimHeatmap(cd8_df, dims = 1:15, cells = 500, balanced = TRUE)

#cd8_df <- JackStraw(cd8_df, num.replicate = 100)
#cd8_df <- ScoreJackStraw(cd8_df, dims = 1:20)
#JackStrawPlot(cd8_df, dims = 1:15)
ElbowPlot(cd8_df, ndims = 30)
# basically decreasing all the time. Maybe 30 would be best though.

#cluster the cells
cd8_df <- FindNeighbors(cd8_df, dims = 1:30)
cd8_df <- FindClusters(cd8_df, resolution = 1.2)
cd8_df <- RunUMAP(cd8_df, dims = 1:30)
DimPlot(cd8_df, reduction = "umap",label=T)
ggsave("umap_clusters_res12.pdf")


#DE markers
cd8_df.markers <- FindAllMarkers(cd8_df, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(cd8_df.markers,file="/scratch/project_2005392/JoonaDawitAnalysisResults/JS/cd8/cd8_df_clutseringMarkers.csv")

maturaatio <- c("CD3D","CD8A", "CD8B","CD4","CCR7","TCF7", "SELL", "S100A4", "GZMA", "GZMB", "GZMM", "GZMH", "GZMK", "PRF1")
FeaturePlot(cd8_df, features = maturaatio)
ggsave("cd8_featureplot_maturation_subclusters.pdf", width = 8, height = 8)
VlnPlot(cd8_df, features = maturaatio)
ggsave("cd8_vlnreplot_maturation_subclusters.pdf", width = 8, height = 8)

exhaustio <- c("TIGIT", "HAVCR2", "PDCD1", "CD244", "LAG3", "TCF7", "CTLA4") #HAVCR2 = 2B4, HAVCR2 = TIM3
FeaturePlot(cd8_df, features = exhaustio)
ggsave("cd8_featureplot_exhaustion_subclusters.pdf", width = 8, height = 8)

VlnPlot(cd8_df, features = exhaustio)
ggsave("cd8_vlnreplot_exhaustion_subclusters.pdf", width = 8, height = 8)



#plots
VlnPlot(cd8_df, features = c("S.Score"))
#no evident differencies in S score
VlnPlot(cd8_df, features = c("G2M.Score"))
#no evident differencies in G2M score

cd8_df.markers %>%
        group_by(cluster) %>%
        slice_max(n = 5, order_by = avg_log2FC) %>% print(n = 100)


#so clusters 12 and 13 are of interest

cd8_df.markers %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top20
DoHeatmap(cd8_df, features = top20$gene) + NoLegend()
ggsave("cd8_heatmap_clusters_top20.pdf", width = 8, height = 8)



# there is a sample-specific batch effect in subclustered CD8s
# Less batch effect on the memory CD8 side, but the naive CD8s appear to be clustered by sample

DimPlot(cd8_df, reduction = "umap",label=F,group.by=c("seurat_clusters","sampleName"),cols=DiscretePalette(15, palette = "glasbey"))
ggsave("cd8_umap_clusters_sampleSpecificBatches.pdf", height = 100, width = 200, units = "mm")

# Checked if there is sample specific batch effect in the B cell reclustering.
# In the B cells, there is better mixing of clusters by different samle names, so is better, and no need to worry about batch effect here
breclust = readRDS("data/msFNAmerged.goodSamples.membcells")
DimPlot(breclust, reduction = "umap",label=F,group.by=c("memnaive_type","sampleName"),cols=DiscretePalette(15, palette = "glasbey"))
ggsave("Umap_B_reclustering_No_sampleSpecificBatches.pdf",height = 100, width = 200, units = "mm")



#### Subclustering again after integrating the CD8 data by donors ####

# since there appears to be a sample-specific batch effect, we decided to do the subclustering after integrating the CD8 data 

cd8_df.List.forIntegration <- SplitObject(cd8_df, split.by = "sampleName")

# normalize and identify variable features for each dataset independently
cd8_df.List.forIntegration <- lapply(X = cd8_df.List.forIntegration, FUN = function(x) {
  x <- NormalizeData(x)
  x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000) 
})

cd8_df.features <- SelectIntegrationFeatures(object.list = cd8_df.List.forIntegration,nfeatures = 2000)

cd8_df.List.forIntegration <- lapply(X = cd8_df.List.forIntegration, FUN = function(x) {
  x <- ScaleData(x, features = cd8_df.features, verbose = FALSE)
  x <- RunPCA(x, features = cd8_df.features, approx=FALSE,verbose = FALSE)
})

cd8_df.selected.anchors <- FindIntegrationAnchors(object.list = cd8_df.List.forIntegration, anchor.features = cd8_df.features)

cd8_df.integrated <- IntegrateData(anchorset = cd8_df.selected.anchors,k.weight = 50)

# Run the standard workflow for visualization and clustering
cd8_df.integrated <- ScaleData(cd8_df.integrated, verbose = FALSE)
cd8_df.integrated <- RunPCA(cd8_df.integrated, npcs = 50, verbose = FALSE)

cd8_df.integrated <- RunUMAP(cd8_df.integrated, reduction = "pca", dims = 1:50)
cd8_df.integrated <- FindNeighbors(cd8_df.integrated, reduction = "pca", dims = 1:50)
cd8_df.integrated <- FindClusters(cd8_df.integrated)

DimPlot(cd8_df.integrated,label=T,cols=DiscretePalette(35, palette = "glasbey"))

DimPlot(cd8_df.integrated,group.by = c("seurat_clusters","sampleName"),cols=DiscretePalette(25, palette = "glasbey"))
ggsave("Umap_cd8_df.integrated_bySampleName.pdf", height = 100, width = 200, units = "mm")

DimPlot(cd8_df.integrated,group.by = c("seurat_clusters","rnaclustersvsAdtManAnn2MaxAnn2"),cols=DiscretePalette(25, palette = "glasbey"))
ggsave("Umap_cd8_df.integrated_byOurAnnotation.pdf", height = 100, width = 200, units = "mm")

DimPlot(cd8_df.integrated,group.by = c("seurat_clusters","immFine"),cols=DiscretePalette(25, palette = "glasbey"))
ggsave("Umap_cd8_df.integrated_byImmFine.pdf", height = 100, width = 200, units = "mm")


# Run the standard workflow for visualization and clustering with higher resolution of 1.4

cd8_df.js <- IntegrateData(anchorset = cd8_df.selected.anchors,k.weight = 50)

cd8_df.js <- ScaleData(cd8_df.js, verbose = FALSE)
cd8_df.js <- RunPCA(cd8_df.js, npcs = 50, verbose = FALSE)

cd8_df.js <- RunUMAP(cd8_df.js, reduction = "pca", dims = 1:50)
cd8_df.js <- FindNeighbors(cd8_df.js, reduction = "pca", dims = 1:50)
cd8_df.js <- FindClusters(cd8_df.js, resolution = 1.4)


#js
DimPlot(cd8_df.js,group.by = c("seurat_clusters","sampleName"),cols=DiscretePalette(25, palette = "glasbey"))
ggsave("Umap_cd8_df.js_bySampleName.pdf", height = 100, width = 200, units = "mm")

DimPlot(cd8_df.js,group.by = c("seurat_clusters","rnaclustersvsAdtManAnn2MaxAnn2"),cols=DiscretePalette(25, palette = "glasbey"))
ggsave("Umap_cd8_df.js_byOurAnnotation.pdf", height = 100, width = 200, units = "mm")

DimPlot(cd8_df.js,group.by = c("seurat_clusters","immFine"),cols=DiscretePalette(25, palette = "glasbey"))
ggsave("Umap_cd8_df.js_byImmFine.pdf", height = 100, width = 200, units = "mm")



DefaultAssay(cd8_df.integrated) <- "RNA"
DefaultAssay(cd8_df.js) <- "RNA"

#DE markers
cd8_df.integrated.markers <- FindAllMarkers(cd8_df.integrated, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cd8_df.integrated.markers.sig = cd8_df.integrated.markers[cd8_df.integrated.markers$p_val_adj < 0.05,]

write.csv(cd8_df.integrated.markers.sig,file="cd8_df.integrated.markers.sig.csv")

cd8_df.integrated.markers.sig %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>% print(n = 100)

#DE markers js
cd8_df.js.markers <- FindAllMarkers(cd8_df.js, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)

cd8_df.js.markers.sig = cd8_df.js.markers[cd8_df.js.markers$p_val_adj < 0.05,]

write.csv(cd8_df.js.markers.sig,file="cd8_df.js.markers.sig.csv")


cd8_df.js.markers.sig %>%
  group_by(cluster) %>%
  slice_max(n = 5, order_by = avg_log2FC) %>% print(n = 100)




#DA tables
cd8VersusSampleCount = table(cd8_df.integrated[["seurat_clusters"]][,1],cd8_df.integrated@meta.data$sampleName)
ctrlSamples <- grepl("CTRL",colnames(cd8VersusSampleCount))
msSamples <- !grepl("CTRL",colnames(cd8VersusSampleCount))

# normalize all samples to a total of 1k cells.
cd8VersusSampleCountNorm <- apply(cd8VersusSampleCount + 1,2, function(x) (x/sum(x)))

# perform unpaired t.test
diffCellAbundance = t(apply(cd8VersusSampleCountNorm,1,function(x) c(t.test(x[msSamples],x[ctrlSamples])$p.value,
                                                                     log2(mean(x[msSamples])/mean(x[ctrlSamples])))))
colnames(diffCellAbundance) <- c("unpaired two-samples t.test.p.value","log2FC(meanMS/meanCTRL)")

fullDiffCellData = cbind(cd8VersusSampleCount,cd8VersusSampleCountNorm,diffCellAbundance)
fullDiffCellData <- fullDiffCellData[order(fullDiffCellData[,which(grepl("p.value",colnames(fullDiffCellData)))],decreasing=F),]

write.csv(fullDiffCellData,"cd8_df.integrated_subtypes_differentialCellTypeAbudance.csv")

#saveRDS
#saveRDS(cd8_df.integrated, file = "/scratch/project_2005392/JoonaDawitAnalysisResults/JS/cd8/cd8_df.integrated")
#cd8_df.integrated <- readRDS(file = "/scratch/project_2005392/JoonaDawitAnalysisResults/JS/cd8/cd8_df.integrated")

#DA tables high resolution clustering by js
cd8VersusSampleCount.js = table(cd8_df.js[["seurat_clusters"]][,1],cd8_df.js@meta.data$sampleName)
ctrlSamples.js <- grepl("CTRL",colnames(cd8VersusSampleCount.js))
msSamples.js <- !grepl("CTRL",colnames(cd8VersusSampleCount.js))
# normalize all samples to a total of 1k cells.
cd8VersusSampleCountNorm.js <- apply(cd8VersusSampleCount.js + 1,2, function(x) (x/sum(x)))
# perform unpaired t.test
diffCellAbundance.js = t(apply(cd8VersusSampleCountNorm.js,1,function(x) c(t.test(x[msSamples.js],x[ctrlSamples.js])$p.value,
                                                                     log2(mean(x[msSamples.js])/mean(x[ctrlSamples.js])))))
colnames(diffCellAbundance.js) <- c("unpaired two-samples t.test.p.value","log2FC(meanMS/meanCTRL)")
fullDiffCellData.js = cbind(cd8VersusSampleCount.js,cd8VersusSampleCountNorm.js,diffCellAbundance.js)
fullDiffCellData.js <- fullDiffCellData.js[order(fullDiffCellData.js[,which(grepl("p.value",colnames(fullDiffCellData.js)))],decreasing=F),]
write.csv(fullDiffCellData.js,"cd8_df.js_subtypes_differentialCellTypeAbudance.csv")

#saveRDS js
#saveRDS(cd8_df.js, file = "/scratch/project_2005392/JoonaDawitAnalysisResults/JS/cd8/cd8_df.js")
#cd8_df.js <- readRDS(file = "/scratch/project_2005392/JoonaDawitAnalysisResults/JS/cd8/cd8_df.js")



# marker plots

maturaatio <- c("CD3D","CD8A", "CD8B","CD4","CCR7","TCF7", "SELL", "S100A4", "GZMA", "GZMB", "GZMM", "GZMH", "GZMK", "PRF1")
FeaturePlot(cd8_df.integrated, features = maturaatio)
ggsave("cd8_df.integrated_featureplot_maturation_subclusters.pdf", width = 8, height = 8)

VlnPlot(cd8_df.integrated, features = maturaatio)
ggsave("cd8_df.integrated_vlnreplot_maturation_subclusters.pdf", width = 8, height = 8)

exhaustio <- c("TIGIT", "HAVCR2", "PDCD1", "CD244", "LAG3", "TCF7", "CTLA4") #HAVCR2 = 2B4, HAVCR2 = TIM3
FeaturePlot(cd8_df.integrated, features = exhaustio)
ggsave("cd8_df.integrated_featureplot_exhaustion_subclusters.pdf", width = 8, height = 8)

VlnPlot(cd8_df.integrated, features = exhaustio)
ggsave("cd8_df.integrated_vlnreplot_exhaustion_subclusters.pdf", width = 8, height = 8)


DefaultAssay(cd8_df.integrated) <- "integrated"

cd8_df.integrated.markers.sig %>%
  group_by(cluster) %>%
  top_n(n = 20, wt = avg_log2FC) -> top10

DoHeatmap(subset(cd8_df.integrated, downsample = 100), features = top10$gene) + NoLegend()
ggsave("cd8_df.integrated_heatmap_clusters_top10.pdf", width = 8, height = 20)

#marker plots high resolution clustering data
FeaturePlot(cd8_df.js, features = maturaatio)
ggsave("cd8_df.js_featureplot_maturation_subclusters.pdf", width = 8, height = 8)
VlnPlot(cd8_df.js, features = maturaatio)
ggsave("cd8_df.js_vlnreplot_maturation_subclusters.pdf", width = 8, height = 8)
FeaturePlot(cd8_df.js, features = exhaustio)
ggsave("cd8_df.js_featureplot_exhaustion_subclusters.pdf", width = 8, height = 8)
VlnPlot(cd8_df.js, features = exhaustio)
ggsave("cd8_df.js_vlnreplot_exhaustion_subclusters.pdf", width = 8, height = 8)

DefaultAssay(cd8_df.js) <- "integrated"
DoHeatmap(subset(cd8_df.js, downsample = 100), features = top10$gene) + NoLegend()
ggsave("cd8_df.js_heatmap_clusters_top10.pdf", width = 8, height = 20)


# Exhausted cluster by cell scoring using the exhaustion markers
# With the markers from the reviewer, cluster 2 and 8 of CD8s have highest exhaustion score, both are memory CD8s
# with our markers, which include TCF7, cluster 6 has highest score, but it is a naive CD8 cluster, likely TCF7 is a marker of naivety as well
# it is probably better to use the markers from the reviewer

exhaustio <- c("TIGIT", "HAVCR2", "PDCD1", "CD244", "LAG3", "TCF7", "CTLA4") #CD244 = 2B4, HAVCR2 = TIM3
ExhaustionMarkersFromReviewer = c("HAVCR2", "PDCD1", "CD244", "LAG3")

Exhaustion_features <- list(EXmarkers1=exhaustio,EXmarkers2=ExhaustionMarkersFromReviewer)


cd8_df.integrated <-  AddModuleScore(object = cd8_df.integrated,assay = "RNA",
                                    features = Exhaustion_features,ctrl = 100,
                                    name = 'Exhaustion_score')

nExhaustionPerCD8Cluster <- cd8_df.integrated@meta.data %>% group_by(seurat_clusters) %>% 
  summarise(n=n(),nEX = sum(Exhaustion_score1 > 0)) %>% mutate(pEX = (nEX/n) * 100) %>% arrange(desc(pEX))

write.csv(nExhaustionPerCD8Cluster,file="nExhaustionPerCD8Cluster.csv")

nExhaustionPerCD8Cluster_ReviewerMarkers <- cd8_df.integrated@meta.data %>% group_by(seurat_clusters) %>% 
  summarise(n=n(),nEX = sum(Exhaustion_score2 > 0)) %>% mutate(pEX = (nEX/n) * 100) %>% arrange(desc(pEX))

write.csv(nExhaustionPerCD8Cluster_ReviewerMarkers,file="nExhaustionPerCD8Cluster_ReviewerMarkers.csv")


FeaturePlot(cd8_df.integrated, label=T,features = "Exhaustion_score2",cols=c("lightblue","blue","Darkblue"))
ggsave("cd8_df.integrated_FeaturePlot_exhaustionMarker2_score.pdf", width = 8, height = 8)

FeaturePlot(cd8_df.integrated, label=T,features = "Exhaustion_score1",cols=c("lightblue","blue","Darkblue"))
ggsave("cd8_df.integrated_FeaturePlot_exhaustionMarker1_score.pdf", width = 8, height = 8)

VlnPlot(subset(cd8_df.integrated, subset = Exhaustion_score2 > 0) , features = "Exhaustion_score2")
ggsave("cd8_df.integrated_vlnreplot_exhaustionMarker2_score.pdf", width = 8, height = 8)


# Compare exhaustion levels between MS and CTRL in CD8s

pCD8_exhaustionLevels = c(0,0,0,0,0,0,0,0,0)
names(pCD8_exhaustionLevels) <- names(table(cd8_df.integrated$sampleName))

cd8_df.integrated.highExhaustionScoreCells = subset(cd8_df.integrated, subset = Exhaustion_score2 > 0)

pEx = table(cd8_df.integrated.highExhaustionScoreCells$sampleName)/table(cd8_df.integrated$sampleName)

pCD8_exhaustionLevels[names(pEx)] <- pEx

pCD8_exhaustionLevels <- pCD8_exhaustionLevels * 100
pCD8_exhaustionLevels <- as.data.frame(pCD8_exhaustionLevels)
pCD8_exhaustionLevels$sampleGroup = sapply(rownames(pCD8_exhaustionLevels),function(jj) ifelse(grepl("CTRL",jj),"CTRL","MS"))

ggplot(pCD8_exhaustionLevels, aes(x=sampleGroup, y=pCD8_exhaustionLevels,fill=sampleGroup)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:2])+
  geom_label(
    label=rownames(pCD8_exhaustionLevels),label.padding = unit(0.08, "lines"), 
    nudge_x = 0.30,nudge_y = -0.05,
    check_overlap = T,size=1
  ) +
  stat_compare_means(aes(group = sampleGroup), label = "p.format", method = "t.test") + 
  theme_classic() + theme(legend.position = "none") +
  ylab("% Exhausted CD8s")

ggsave("cd8_df.integrated_percentOfExhaustedCD8s_exhaustionMarkers2Score_betweenConditions.pdf", width = 2, height = 4)


## Check just within cluster 2, which is the exhausted cluster

pCD8_exhaustionLevels_cls2 = c(0,0,0,0,0,0,0,0,0)
names(pCD8_exhaustionLevels_cls2) <- names(table(cd8_df.integrated$sampleName))

cd8_df.integrated.highExhaustionScoreCells.cls2 = subset(cd8_df.integrated, subset = Exhaustion_score2 > 0 & seurat_clusters==2)

cd8_df.integrated.cls2 = subset(cd8_df.integrated, subset = seurat_clusters==2)


pEx = table(cd8_df.integrated.highExhaustionScoreCells.cls2$sampleName)/table(cd8_df.integrated.cls2$sampleName)

pCD8_exhaustionLevels_cls2[names(pEx)] <- pEx

pCD8_exhaustionLevels_cls2 <- pCD8_exhaustionLevels_cls2 * 100
pCD8_exhaustionLevels_cls2 <- as.data.frame(pCD8_exhaustionLevels_cls2)
pCD8_exhaustionLevels_cls2$sampleGroup = sapply(rownames(pCD8_exhaustionLevels_cls2),function(jj) ifelse(grepl("CTRL",jj),"CTRL","MS"))

ggplot(pCD8_exhaustionLevels_cls2, aes(x=sampleGroup, y=pCD8_exhaustionLevels_cls2,fill=sampleGroup)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:2])+
  geom_label(
    label=rownames(pCD8_exhaustionLevels_cls2),label.padding = unit(0.08, "lines"), 
    nudge_x = 0.30,nudge_y = -0.05,
    check_overlap = T,size=1
  ) +
  stat_compare_means(aes(group = sampleGroup), label = "p.format", method = "t.test") + 
  theme_classic() + theme(legend.position = "none") +
  ylab("% Exhausted in cluster 2 of CD8s")

ggsave("cd8_df.integrated_percentOfExhaustedCD8sInCluster2_exhaustionMarkers2Score_betweenConditions.pdf", width = 2, height = 4)



## Evaluation of exhaustion after the higher resolution clustering 

# Joona suggested cluster 9 to be a candidate for exhausted CD8 cluster. Here we check with the scoring approach.

# Exhausted cluster by cell scoring using the exhaustion markers
# With the markers from the reviewer, cluster 9 of CD8s has highest exhaustion score (much higher than the others). But with the 
# our list of markers for exhaustion (exhaustio below), cluster 13 has the highest score, this is definitely related to the fact that
# cluster 13 expresses TCF7 higher than other clusters just like cluster 8 and 11, all of which are scored higher than cluster 9 with this list
# of markers. TCF7 is also a marker of naivety or stemness. are these some how related ?  
# In any case, We decided to use the markers from the reviewer as it clearly indicates cluster 9 is the exhausted cluster.

exhaustio <- c("TIGIT", "HAVCR2", "PDCD1", "CD244", "LAG3", "TCF7", "CTLA4") #CD244 = 2B4, HAVCR2 = TIM3
ExhaustionMarkersFromReviewer = c("HAVCR2", "PDCD1", "CD244", "LAG3")

Exhaustion_features <- list(EXmarkers1=exhaustio,EXmarkers2=ExhaustionMarkersFromReviewer)


cd8_df.js <-  AddModuleScore(object = cd8_df.js,assay = "RNA",
                                     features = Exhaustion_features,ctrl = 100,
                                     name = 'Exhaustion_score')

nExhaustionPerCD8Cluster.js <- cd8_df.js@meta.data %>% group_by(seurat_clusters) %>% 
  summarise(n=n(),nEX = sum(Exhaustion_score1 > 0)) %>% mutate(pEX = (nEX/n) * 100) %>% arrange(desc(pEX))

write.csv(nExhaustionPerCD8Cluster.js,file="nExhaustionPerCD8Cluster.js.csv")

nExhaustionPerCD8Cluster_ReviewerMarkers.js <- cd8_df.js@meta.data %>% group_by(seurat_clusters) %>% 
  summarise(n=n(),nEX = sum(Exhaustion_score2 > 0)) %>% mutate(pEX = (nEX/n) * 100) %>% arrange(desc(pEX))

write.csv(nExhaustionPerCD8Cluster_ReviewerMarkers.js,file="nExhaustionPerCD8Cluster_ReviewerMarkers.js.csv")


FeaturePlot(cd8_df.js, label=T,features = "Exhaustion_score2",cols=c("lightblue","blue","Darkblue"))
ggsave("cd8_df.js_FeaturePlot_exhaustionMarker2_score.pdf", width = 8, height = 8)

FeaturePlot(cd8_df.js, label=T,features = "Exhaustion_score1",cols=c("lightblue","blue","Darkblue"))
ggsave("cd8_df.js.integrated_FeaturePlot_exhaustionMarker1_score.pdf", width = 8, height = 8)

VlnPlot(subset(cd8_df.js, subset = Exhaustion_score2 > 0) , features = "Exhaustion_score2")
ggsave("cd8_df.js_vlnreplot_exhaustionMarker2_score.pdf", width = 8, height = 8)


DoHeatmap(subset(cd8_df.js, downsample = 100), features = exhaustio) 
ggsave("cd8_df.js_heatmap_clusters_exhaustionmarkers.pdf", width = 10, height = 12)





# Compare exhaustion levels between MS and CTRL in CD8s

pCD8_exhaustionLevels.js = c(0,0,0,0,0,0,0,0,0)
names(pCD8_exhaustionLevels.js) <- names(table(cd8_df.js$sampleName))

cd8_df.js.highExhaustionScoreCells = subset(cd8_df.js, subset = Exhaustion_score2 > 0)

pEx.js = table(cd8_df.js.highExhaustionScoreCells$sampleName)/table(cd8_df.js$sampleName)

pCD8_exhaustionLevels.js[names(pEx.js)] <- pEx.js

pCD8_exhaustionLevels.js <- pCD8_exhaustionLevels.js * 100
pCD8_exhaustionLevels.js <- as.data.frame(pCD8_exhaustionLevels.js)
pCD8_exhaustionLevels.js$sampleGroup = sapply(rownames(pCD8_exhaustionLevels.js),function(jj) ifelse(grepl("CTRL",jj),"CTRL","MS"))

ggplot(pCD8_exhaustionLevels.js, aes(x=sampleGroup, y=pCD8_exhaustionLevels.js,fill=sampleGroup)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:2])+
  geom_label(
    label=rownames(pCD8_exhaustionLevels.js),label.padding = unit(0.08, "lines"), 
    nudge_x = 0.30,nudge_y = -0.05,
    check_overlap = T,size=1
  ) +
  stat_compare_means(aes(group = sampleGroup), label = "p.format", method = "t.test") + 
  theme_classic() + theme(legend.position = "none") +
  ylab("% Exhausted CD8s")

ggsave("cd8_df.js_percentOfExhaustedCD8s_exhaustionMarkers2Score_betweenConditions.pdf", width = 2, height = 4)


## Check just within cluster 9, which is the exhausted cluster

pCD8_exhaustionLevels_cls9.js = c(0,0,0,0,0,0,0,0,0)
names(pCD8_exhaustionLevels_cls9.js) <- names(table(cd8_df.js$sampleName))

cd8_df.js.highExhaustionScoreCells.cls9 = subset(cd8_df.js, subset = Exhaustion_score2 > 0 & seurat_clusters==9)

cd8_df.js.cls9 = subset(cd8_df.js, subset = seurat_clusters==9)


pEx.js9 = table(cd8_df.js.highExhaustionScoreCells.cls9$sampleName)/table(cd8_df.js.cls9$sampleName)

pCD8_exhaustionLevels_cls9.js[names(pEx.js9)] <- pEx.js9

pCD8_exhaustionLevels_cls9.js <- pCD8_exhaustionLevels_cls9.js * 100
pCD8_exhaustionLevels_cls9.js <- as.data.frame(pCD8_exhaustionLevels_cls9.js)
pCD8_exhaustionLevels_cls9.js$sampleGroup = sapply(rownames(pCD8_exhaustionLevels_cls9.js),function(jj) ifelse(grepl("CTRL",jj),"CTRL","MS"))

ggplot(pCD8_exhaustionLevels_cls9.js, aes(x=sampleGroup, y=pCD8_exhaustionLevels_cls9.js,fill=sampleGroup)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:2])+
  geom_label(
    label=rownames(pCD8_exhaustionLevels_cls9.js),label.padding = unit(0.08, "lines"), 
    nudge_x = 0.30,nudge_y = -0.05,
    check_overlap = T,size=1
  ) +
  stat_compare_means(aes(group = sampleGroup), label = "p.format", method = "t.test") + 
  theme_classic() + theme(legend.position = "none") +
  ylab("% Exhausted in cluster 9 of CD8s")

ggsave("cd8_df.js_percentOfExhaustedCD8sInCluster9_exhaustionMarkers2Score_betweenConditions.pdf", width = 2, height = 4)


## % of exhausted cells as: number of cluster 9 cells in sample / total number of CD8 cells in sample

sizeOfCluster9_insamples = table(cd8_df.js$sampleName,cd8_df.js$seurat_clusters)[,10]
nCd8_cells_insamples = table(cd8_df.js$sampleName)

sizeOfCluster9_insamples/nCd8_cells_insamples

compareSizeOfCD8subclusters <- function(srtobjct,clusterName){
  
  
  sizeOfCluster_insamples = table(srtobjct@meta.data[srtobjct@meta.data$seurat_clusters==clusterName,]$sampleName)
  nCd8_cells_insamples = table(srtobjct$sampleName)
  
  propOfClustInSample = 100 * sizeOfCluster_insamples/nCd8_cells_insamples
  

  propOfClustInSample.df <- as.data.frame(propOfClustInSample)
  rownames(propOfClustInSample.df) <- propOfClustInSample.df[,1]
  propOfClustInSample.df <- propOfClustInSample.df[,-1,drop=F]
  
  propOfClustInSample.df$sampleGroup = sapply(rownames(propOfClustInSample.df),function(jj) ifelse(grepl("CTRL",jj),"CTRL","MS"))
  
  ggplot(propOfClustInSample.df, aes(x=sampleGroup, y=Freq,fill=sampleGroup)) + 
    geom_boxplot(outlier.colour = NA) +
    geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:2])+
    geom_label(
      label=rownames(propOfClustInSample.df),label.padding = unit(0.08, "lines"),size=2
    ) +
    stat_compare_means(aes(group = sampleGroup), label = "p.format", method = "t.test") + 
    theme_classic() + theme(legend.position = "none") +
    ylab(paste("Proportion of subcluster",clusterName,"(% of all CD8 cells)"))
  
  ggsave(paste0("/scratch/project_2005392/JoonaDawitAnalysisResults/JS/cd8/cd8_df.js_percentOf_subcluster_",clusterName,"_betweenConditions.pdf"), width = 2, height = 4)
  
  
  
  
}

for(cd8cls in levels(cd8_df.js$seurat_clusters)){
  
  compareSizeOfCD8subclusters(srtobjct=cd8_df.js,clusterName=cd8cls)
}


## TCR sharing between the CD8 subclusters ####

# where do the clones in the "exhausted CD8 cluster" arise from most ?

cellsWithTCRs = rownames(cd8_df.js@meta.data)[!is.na(cd8_df.js@meta.data$CTstrict)]

cd8_df.js.withTCRs = cd8_df.js[,cellsWithTCRs]

clusterClonotypes_vs_sharedclusters = NULL

cd8_cls = levels(Idents(cd8_df.js.withTCRs))

for(sm in unique(cd8_df.js.withTCRs$sampleName)){
  
  cd8_df.js.withTCRs.sm = subset(cd8_df.js.withTCRs,subset = sampleName == sm)
  
  for(clsid in cd8_cls){
    cd8_df.js.withTCRs.cls1 = cd8_df.js.withTCRs.sm@meta.data[cd8_df.js.withTCRs.sm@meta.data$seurat_clusters==clsid,]
    
    for(clsid2 in cd8_cls){
      cd8_df.js.withTCRs.cls2 = cd8_df.js.withTCRs.sm@meta.data[cd8_df.js.withTCRs.sm@meta.data$seurat_clusters==clsid2,]
      
      clsid1clonesIn2 = sapply(cd8_df.js.withTCRs.cls1$CTstrict, function(x) sum(cd8_df.js.withTCRs.cls2$CTstrict == x))
      numberOfOverlaps = sum(clsid1clonesIn2 > 0)
      maxSizeOfCellsWithSharedClonotype = max(as.numeric(names(table(clsid1clonesIn2))))
      
      overlapRes = c(sm,clsid,clsid2,nrow(cd8_df.js.withTCRs.cls1),length(table(cd8_df.js.withTCRs.cls1$CTstrict)),
                     nrow(cd8_df.js.withTCRs.cls2),length(table(cd8_df.js.withTCRs.cls2$CTstrict)),
                     numberOfOverlaps,maxSizeOfCellsWithSharedClonotype)
      names(overlapRes) = c("sampleName","CD8cls1","CD8cls2","CD8cls1_ncells","CD8cls1_nClonotypes",
                            "CD8cls2_ncells","CD8cls2_nClonotypes","cls1_cls2_nCellsWithSharedabTCRs","highestDetectednCellsWithSameabTCRs")
      
      clusterClonotypes_vs_sharedclusters = rbind(clusterClonotypes_vs_sharedclusters,overlapRes)
      
    }
    
    
  }
  
}

clusterClonotypes_vs_sharedclusters = as.data.frame(clusterClonotypes_vs_sharedclusters)

# only cluster 9 sharings 
clusterClonotypes_vs_sharedclusters_cls9 = clusterClonotypes_vs_sharedclusters[clusterClonotypes_vs_sharedclusters$CD8cls1==9,]

# very clearly, clonotypes in the exhausted cluster appears to originate largely from cluster 6 both in MS and CTRLs, and more so in MS patient
# As cluster 6 is the most increased in MS patients, it is clearly the source of the exhausted cells in MS,
# but also in controls since it is an active population, probably the most active.

# using heatmap to show this result

ClsCombinations = unique(paste(clusterClonotypes_vs_sharedclusters_cls9$CD8cls1,clusterClonotypes_vs_sharedclusters_cls9$CD8cls2,sep="_"))
smnames = unique(clusterClonotypes_vs_sharedclusters_cls9$sampleName)

clusterClonotypes_vs_sharedclusters_cls9_forheatmap = NULL
clusterClonotypes_vs_sharedclusters_cls9_forheatmap2 = NULL # highest size of shared clonotype


for(sm in smnames){
  
  clusterClonotypes_vs_sharedclusters_cls9_sm = clusterClonotypes_vs_sharedclusters_cls9[clusterClonotypes_vs_sharedclusters_cls9$sampleName == sm,]
  
  ovNames = paste(clusterClonotypes_vs_sharedclusters_cls9_sm$CD8cls1,clusterClonotypes_vs_sharedclusters_cls9_sm$CD8cls2,sep="_")
  
  clusterClonotypes_vs_sharedclusters_cls9_sm$nsharedPercOfCls9 = 100 * as.numeric(clusterClonotypes_vs_sharedclusters_cls9_sm$cls1_cls2_nCellsWithSharedabTCRs)/as.numeric(clusterClonotypes_vs_sharedclusters_cls9_sm$CD8cls1_ncells)
  
  clusterClonotypes_vs_sharedclusters_cls9_sm$nsharedPercOfBothcls = 100 * 2 * (as.numeric(clusterClonotypes_vs_sharedclusters_cls9_sm$cls1_cls2_nCellsWithSharedabTCRs)/(as.numeric(clusterClonotypes_vs_sharedclusters_cls9_sm$CD8cls1_ncells) + as.numeric(clusterClonotypes_vs_sharedclusters_cls9_sm$CD8cls2_ncells)))
  
  
  ovbetweenCls = clusterClonotypes_vs_sharedclusters_cls9_sm$nsharedPercOfCls9
  
  ovhighestSharedCls = clusterClonotypes_vs_sharedclusters_cls9_sm$highestDetectednCellsWithSameabTCRs
  
  
  names(ovbetweenCls) = ovNames
  
  names(ovhighestSharedCls) = ovNames
  
  clusterClonotypes_vs_sharedclusters_cls9_forheatmap = rbind(clusterClonotypes_vs_sharedclusters_cls9_forheatmap,
                                                              c(as.numeric(clusterClonotypes_vs_sharedclusters_cls9_sm$CD8cls1_ncells)[1],as.numeric(ovbetweenCls[ClsCombinations])))
  
  
  clusterClonotypes_vs_sharedclusters_cls9_forheatmap2 = rbind(clusterClonotypes_vs_sharedclusters_cls9_forheatmap2,
                                                               c(as.numeric(clusterClonotypes_vs_sharedclusters_cls9_sm$CD8cls1_ncells)[1],as.numeric(ovhighestSharedCls[ClsCombinations])))
  
}

rownames(clusterClonotypes_vs_sharedclusters_cls9_forheatmap) <- smnames
colnames(clusterClonotypes_vs_sharedclusters_cls9_forheatmap) <- c("cluster9_ncells",ClsCombinations)

rownames(clusterClonotypes_vs_sharedclusters_cls9_forheatmap2) <- smnames
colnames(clusterClonotypes_vs_sharedclusters_cls9_forheatmap2) <- c("cluster9_ncells",ClsCombinations)



pheatmap::pheatmap(clusterClonotypes_vs_sharedclusters_cls9_forheatmap,scale = "none")






#### clonality in CD8 subclusters ####

DefaultAssay(cd8_df.integrated) <- "RNA"

tcrCells = rownames(msFNA.merged.tcr@meta.data[!is.na(msFNA.merged.tcr@meta.data$barcode),])

tcrCells_cd8s = colnames(cd8_df.integrated)[colnames(cd8_df.integrated) %in% tcrCells]

msFNA.merged.CD8.tcrData = subset(msFNA.merged.tcr, subset = barcode %in% tcrCells_cd8s)

cd8_cellTypesWithTCRs <- unique(cd8_df.integrated@meta.data$seurat_clusters)
ResultDir = "cd8/"

for(ctcr in cd8_cellTypesWithTCRs){
  
  print(ctcr)
  
  cd8_df_ctcr = cd8_df.integrated@meta.data[cd8_df.integrated@meta.data$seurat_clusters==ctcr,]
  
  msFNA.merged.CD8.tcrTempData = subset(msFNA.merged.CD8.tcrData,subset = barcode %in% rownames(cd8_df_ctcr))
  combined_MS5_tcrcontigList_tcrTempData <- expression2List(msFNA.merged.CD8.tcrTempData,group="sampleName")
  
  write.to = paste0(ResultDir,ctcr,"_CD8_subcluster_TCRclonality/")
  dir.create(write.to)

  percOfuniqd = quantContig(combined_MS5_tcrcontigList_tcrTempData, cloneCall="CTstrict", scale = T,exportTable = T)
  percOfuniqd.ttest = tryCatch({t.test(percOfuniqd$scaled[grepl("MS",percOfuniqd$values)],percOfuniqd$scaled[grepl("CTRL",percOfuniqd$values)])},
                               error=function(e) {
                                 return(NA)
                               })
  
  write.csv(percOfuniqd,paste0(write.to,ctcr,"_percentOfUniqueClonotypes.csv"))
  if(class(percOfuniqd.ttest) == "htest")
    write.csv(percOfuniqd.ttest$p.value,paste0(write.to,ctcr,"_percentOfUniqueClonotypes_ttestPval.csv"))
  
  percOfuniqd$sampleGroup <- ifelse(grepl("MS",percOfuniqd$values),"MS","CTRL")
  percOfuniqd$sampleGroup <- as.factor(percOfuniqd$sampleGroup)
  
  
  percOfuniqd_statp  = tryCatch({compare_means(scaled ~ sampleGroup, percOfuniqd,method="t.test",p.adjust.method="BH")},
                                error=function(e) {
                                  return(NA)
                                })
  
  if(class(percOfuniqd_statp)[1] == "tbl_df"){
    stat.test <- percOfuniqd_statp %>% mutate(y.position = rep(max(percOfuniqd$scaled) + 0.5,1))
    
    
    g <- ggplot(data = percOfuniqd, aes(x = sampleGroup,y = scaled,fill=sampleGroup))  
    g + geom_boxplot(outlier.colour = NA) +  
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") + 
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Percent of Unique Clonotypes (AA)") + ggtitle(ctcr) +  NoLegend() + 
      
      # significance indicated with stars
      #stat_pvalue_manual(stat.test, x = "timePoint", hide.ns = TRUE,label = "p.signif",label.size=6,inherit.aes = FALSE,tip.length = 0)
      
      # significance with pvalues
      stat_pvalue_manual(stat.test, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)
    
    ggsave(paste0(write.to,ctcr,"_PercentageOfUniqueAAClonotypes_ForPaper.pdf"), width = 2, height = 4)
    
  }else{
    
    
    g <- ggplot(data = percOfuniqd, aes(x = sampleGroup,y = scaled,fill=sampleGroup))  
    g + geom_boxplot(outlier.colour = NA) +  
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") + 
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Percent of Unique Clonotypes (AA)") + ggtitle(ctcr) + NoLegend() 
    
    ggsave(paste0(write.to,ctcr,"_PercentageOfUniqueAAClonotypes_ForPaper.pdf"), width = 2, height = 4)
    
  }
  
  
  
  
  
  # evaluate diversity
  # The shannon diversity is not normalized, needs to be normalized perhaps
  
  clonalDiversity(combined_MS5_tcrcontigList_tcrTempData,cloneCall = "CTstrict",group="sampleGroup",exportTable = F)
  #ggsave(paste0(write.to,ctcr,"_tcrdiversity.pdf"), height = 210, width = 320, units = "mm")
  ggsave(paste0(write.to,ctcr,"_tcrdiversity.pdf"), width = 2, height = 4)
  
  diversitytcrClones = clonalDiversity(combined_MS5_tcrcontigList_tcrTempData,cloneCall = "CTstrict",group="sampleGroup",exportTable = T)
  diversitytcrClones$sample=rownames(diversitytcrClones)
  
  diversitytcrClones.ttest = tryCatch({t.test(diversitytcrClones$Shannon[grepl("MS",diversitytcrClones$sample)],diversitytcrClones$Shannon[grepl("CTRL",diversitytcrClones$sample)])},
                                      error=function(e) {
                                        return(NA)
                                      })
  
  
  #ggplot(diversitytcrClones,aes(x=sampleGroup,y=Shannon,color=sampleGroup)) + 
  # geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + theme_bw() + labs(y="Diversity (Shannon)")
  
  diversitytcrClones$sampleGroup <- as.factor(diversitytcrClones$sampleGroup)
  
  diversitytcrClones_statp  = tryCatch({compare_means(Shannon ~ sampleGroup, diversitytcrClones,method="t.test",p.adjust.method="BH")},
                                       error=function(e) {
                                         return(NA)
                                       })
  
  if(class(diversitytcrClones_statp)[1] == "tbl_df"){
    
    stat.test <- diversitytcrClones_statp %>% mutate(y.position = rep(max(diversitytcrClones$Shannon) + 0.001,1))
    
    g <- ggplot(data = diversitytcrClones, aes(x = sampleGroup,y = Shannon,fill=sampleGroup))  
    g + geom_boxplot(outlier.colour = NA) + 
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") + 
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Shannon diversity") + ggtitle(ctcr) +  NoLegend() + 
      
      # significance indicated with stars
      #stat_pvalue_manual(stat.test, x = "timePoint", hide.ns = TRUE,label = "p.signif",label.size=6,inherit.aes = FALSE,tip.length = 0)
      
      # significance with pvalues
      stat_pvalue_manual(stat.test, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)
    ggsave(paste0(write.to,ctcr,"_tcrdiversity_Shannon.pdf"), width = 2, height = 4)
    
    
  }else{
    
    g <- ggplot(data = diversitytcrClones, aes(x = sampleGroup,y = Shannon,fill=sampleGroup))  
    g + geom_boxplot(outlier.colour = NA) +  
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") + 
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Shannon diversity") + ggtitle(ctcr) + NoLegend() 
    
    ggsave(paste0(write.to,ctcr,"_tcrdiversity_Shannon.pdf"),width = 2, height = 4)
    
  }
  
  
  
  # by downsampling
  
  repsTemp = combined_MS5_tcrcontigList_tcrTempData
  
  downSampleTo=findSmallestNumberOfReads(repsTemp) # smallest repertoire MS002 and has size of 5258 clonotypes
  
  if(downSampleTo < 10)
    next
  
  MSFNA_msTCR_Diversity_downsampledtemp <- estimateDiversity(repsTemp[grepl("MS",names(repsTemp))],downsampleSize=downSampleTo,useSeq="AA",nResamples=100)
  
  MSFNA_ctrlTCR_Diversity_downsampledtemp <- estimateDiversity(repsTemp[grepl("CTRL",names(repsTemp))],downsampleSize=downSampleTo,useSeq="AA",nResamples=100)
  
  # plot diversity dot plot
  diversityValuesAAtemp = data.frame(diversity=c(MSFNA_msTCR_Diversity_downsampledtemp,MSFNA_ctrlTCR_Diversity_downsampledtemp),
                                     sampleGroup=c(rep("MS",length(MSFNA_msTCR_Diversity_downsampledtemp)),
                                                   rep("CTRL",length(MSFNA_ctrlTCR_Diversity_downsampledtemp)))
  )
  
  diversityValuesAAtemp$sampleGroup <- as.factor(diversityValuesAAtemp$sampleGroup)
  
  diversityValuesAAtemp_statp  = compare_means(diversity ~ sampleGroup, diversityValuesAAtemp,method="t.test",p.adjust.method="BH")
  
  stat.test <- diversityValuesAAtemp_statp %>%
    mutate(y.position = rep(max(diversityValuesAAtemp$diversity) + 0.01,1))
  
  
  g <- ggplot(data = diversityValuesAAtemp, aes(x = sampleGroup,y = diversity,fill=sampleGroup))  
  g + geom_boxplot(outlier.colour = NA) +  
    #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
    #scale_color_manual(values=c("lightblue","Brown")) + 
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
    geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
    theme_classic() + theme(legend.position = "none") + 
    theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
    ylab("Normalized Shannon Entropy") +  NoLegend() + 
    
    # significance indicated with stars
    #stat_pvalue_manual(stat.test, x = "timePoint", hide.ns = TRUE,label = "p.signif",label.size=6,inherit.aes = FALSE,tip.length = 0)
    
    # significance with pvalues
    stat_pvalue_manual(stat.test, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)
  
  ggsave(paste0(write.to,ctcr,"_scRepAn_diversityAnalysisAAClonotype_after100Downsampling.pdf"), width = 2, height = 4)
  
  
  
}


## Clonality comparison for the clusters at higher resolution (Joona's clustering with resolution 1.4)
## for seurat object: cd8_df.js


DefaultAssay(cd8_df.js) <- "RNA"

tcrCells = rownames(msFNA.merged.tcr@meta.data[!is.na(msFNA.merged.tcr@meta.data$barcode),])

tcrCells_cd8s = colnames(cd8_df.js)[colnames(cd8_df.js) %in% tcrCells]

msFNA.merged.CD8.tcrData = subset(msFNA.merged.tcr, subset = barcode %in% tcrCells_cd8s)

cd8_cellTypesWithTCRs.js <- unique(cd8_df.js@meta.data$seurat_clusters)
ResultDir = "cd8/"

for(ctcr in cd8_cellTypesWithTCRs.js){
  
  print(ctcr)
  
  cd8_df_ctcr = cd8_df.js@meta.data[cd8_df.js@meta.data$seurat_clusters==ctcr,]
  
  msFNA.merged.CD8.tcrTempData = subset(msFNA.merged.CD8.tcrData,subset = barcode %in% rownames(cd8_df_ctcr))
  combined_MS5_tcrcontigList_tcrTempData <- expression2List(msFNA.merged.CD8.tcrTempData,group="sampleName")
  
  write.to = paste0(ResultDir,ctcr,"_CD8_js_subcluster_TCRclonality/")
  dir.create(write.to)
  
  percOfuniqd = quantContig(combined_MS5_tcrcontigList_tcrTempData, cloneCall="CTstrict", scale = T,exportTable = T)
  percOfuniqd.ttest = tryCatch({t.test(percOfuniqd$scaled[grepl("MS",percOfuniqd$values)],percOfuniqd$scaled[grepl("CTRL",percOfuniqd$values)])},
                               error=function(e) {
                                 return(NA)
                               })
  
  write.csv(percOfuniqd,paste0(write.to,ctcr,"_percentOfUniqueClonotypes.csv"))
  if(class(percOfuniqd.ttest) == "htest")
    write.csv(percOfuniqd.ttest$p.value,paste0(write.to,ctcr,"_percentOfUniqueClonotypes_ttestPval.csv"))
  
  percOfuniqd$sampleGroup <- ifelse(grepl("MS",percOfuniqd$values),"MS","CTRL")
  percOfuniqd$sampleGroup <- as.factor(percOfuniqd$sampleGroup)
  
  
  percOfuniqd_statp  = tryCatch({compare_means(scaled ~ sampleGroup, percOfuniqd,method="t.test",p.adjust.method="BH")},
                                error=function(e) {
                                  return(NA)
                                })
  
  if(class(percOfuniqd_statp)[1] == "tbl_df"){
    stat.test <- percOfuniqd_statp %>% mutate(y.position = rep(max(percOfuniqd$scaled) + 0.5,1))
    
    
    g <- ggplot(data = percOfuniqd, aes(x = sampleGroup,y = scaled,fill=sampleGroup))  
    g + geom_boxplot(outlier.colour = NA) +  
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") + 
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Percent of Unique Clonotypes (AA)") + ggtitle(ctcr) +  NoLegend() + 
      
      # significance indicated with stars
      #stat_pvalue_manual(stat.test, x = "timePoint", hide.ns = TRUE,label = "p.signif",label.size=6,inherit.aes = FALSE,tip.length = 0)
      
      # significance with pvalues
      stat_pvalue_manual(stat.test, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)
    
    ggsave(paste0(write.to,ctcr,"_PercentageOfUniqueAAClonotypes_ForPaper.pdf"), width = 2, height = 4)
    
  }else{
    
    
    g <- ggplot(data = percOfuniqd, aes(x = sampleGroup,y = scaled,fill=sampleGroup))  
    g + geom_boxplot(outlier.colour = NA) +  
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") + 
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Percent of Unique Clonotypes (AA)") + ggtitle(ctcr) + NoLegend() 
    
    ggsave(paste0(write.to,ctcr,"_PercentageOfUniqueAAClonotypes_ForPaper.pdf"), width = 2, height = 4)
    
  }
  
  
  
  
  
  # evaluate diversity
  # The shannon diversity is not normalized, needs to be normalized perhaps
  
  clonalDiversity(combined_MS5_tcrcontigList_tcrTempData,cloneCall = "CTstrict",group="sampleGroup",exportTable = F)
  #ggsave(paste0(write.to,ctcr,"_tcrdiversity.pdf"), height = 210, width = 320, units = "mm")
  ggsave(paste0(write.to,ctcr,"_tcrdiversity.pdf"), width = 2, height = 4)
  
  diversitytcrClones = clonalDiversity(combined_MS5_tcrcontigList_tcrTempData,cloneCall = "CTstrict",group="sampleGroup",exportTable = T)
  diversitytcrClones$sample=rownames(diversitytcrClones)
  
  diversitytcrClones.ttest = tryCatch({t.test(diversitytcrClones$Shannon[grepl("MS",diversitytcrClones$sample)],diversitytcrClones$Shannon[grepl("CTRL",diversitytcrClones$sample)])},
                                      error=function(e) {
                                        return(NA)
                                      })
  
  
  #ggplot(diversitytcrClones,aes(x=sampleGroup,y=Shannon,color=sampleGroup)) + 
  # geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + theme_bw() + labs(y="Diversity (Shannon)")
  
  diversitytcrClones$sampleGroup <- as.factor(diversitytcrClones$sampleGroup)
  
  diversitytcrClones_statp  = tryCatch({compare_means(Shannon ~ sampleGroup, diversitytcrClones,method="t.test",p.adjust.method="BH")},
                                       error=function(e) {
                                         return(NA)
                                       })
  
  if(class(diversitytcrClones_statp)[1] == "tbl_df"){
    
    stat.test <- diversitytcrClones_statp %>% mutate(y.position = rep(max(diversitytcrClones$Shannon) + 0.001,1))
    
    g <- ggplot(data = diversitytcrClones, aes(x = sampleGroup,y = Shannon,fill=sampleGroup))  
    g + geom_boxplot(outlier.colour = NA) + 
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") + 
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Shannon diversity") + ggtitle(ctcr) +  NoLegend() + 
      
      # significance indicated with stars
      #stat_pvalue_manual(stat.test, x = "timePoint", hide.ns = TRUE,label = "p.signif",label.size=6,inherit.aes = FALSE,tip.length = 0)
      
      # significance with pvalues
      stat_pvalue_manual(stat.test, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)
    ggsave(paste0(write.to,ctcr,"_tcrdiversity_Shannon.pdf"), width = 2, height = 4)
    
    
  }else{
    
    g <- ggplot(data = diversitytcrClones, aes(x = sampleGroup,y = Shannon,fill=sampleGroup))  
    g + geom_boxplot(outlier.colour = NA) +  
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") + 
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Shannon diversity") + ggtitle(ctcr) + NoLegend() 
    
    ggsave(paste0(write.to,ctcr,"_tcrdiversity_Shannon.pdf"),width = 2, height = 4)
    
  }
  
  
  
  # by downsampling
  
  repsTemp = combined_MS5_tcrcontigList_tcrTempData
  
  downSampleTo=findSmallestNumberOfReads(repsTemp) # smallest repertoire MS002 and has size of 5258 clonotypes
  
  if(downSampleTo < 10)
    next
  
  MSFNA_msTCR_Diversity_downsampledtemp <- estimateDiversity(repsTemp[grepl("MS",names(repsTemp))],downsampleSize=downSampleTo,useSeq="AA",nResamples=100)
  
  MSFNA_ctrlTCR_Diversity_downsampledtemp <- estimateDiversity(repsTemp[grepl("CTRL",names(repsTemp))],downsampleSize=downSampleTo,useSeq="AA",nResamples=100)
  
  # plot diversity dot plot
  diversityValuesAAtemp = data.frame(diversity=c(MSFNA_msTCR_Diversity_downsampledtemp,MSFNA_ctrlTCR_Diversity_downsampledtemp),
                                     sampleGroup=c(rep("MS",length(MSFNA_msTCR_Diversity_downsampledtemp)),
                                                   rep("CTRL",length(MSFNA_ctrlTCR_Diversity_downsampledtemp)))
  )
  
  diversityValuesAAtemp$sampleGroup <- as.factor(diversityValuesAAtemp$sampleGroup)
  
  diversityValuesAAtemp_statp  = compare_means(diversity ~ sampleGroup, diversityValuesAAtemp,method="t.test",p.adjust.method="BH")
  
  stat.test <- diversityValuesAAtemp_statp %>%
    mutate(y.position = rep(max(diversityValuesAAtemp$diversity) + 0.01,1))
  
  
  g <- ggplot(data = diversityValuesAAtemp, aes(x = sampleGroup,y = diversity,fill=sampleGroup))  
  g + geom_boxplot(outlier.colour = NA) +  
    #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
    #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
    #scale_color_manual(values=c("lightblue","Brown")) + 
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
    geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
    theme_classic() + theme(legend.position = "none") + 
    theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
    ylab("Normalized Shannon Entropy") +  NoLegend() + 
    
    # significance indicated with stars
    #stat_pvalue_manual(stat.test, x = "timePoint", hide.ns = TRUE,label = "p.signif",label.size=6,inherit.aes = FALSE,tip.length = 0)
    
    # significance with pvalues
    stat_pvalue_manual(stat.test, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)
  
  ggsave(paste0(write.to,ctcr,"_scRepAn_diversityAnalysisAAClonotype_after100Downsampling.pdf"), width = 2, height = 4)
  
  
 
}




#### Compare frequency distribution of CD8 subclusters in pooled data ####

# for the high resolution clustering of integrated Cd8 data

cellCountsInMS = table(cd8_df.js@meta.data[cd8_df.js@meta.data$sampleGroup=="MS",]$seurat_clusters)
cellCountsInCTRL = table(cd8_df.js@meta.data[cd8_df.js@meta.data$sampleGroup=="CTRL",]$seurat_clusters)

CD8_subset_cellFrequenyData <- as.table(rbind(cellCountsInMS, cellCountsInCTRL))
dimnames(CD8_subset_cellFrequenyData) <- list(sampleGroup = c("MS", "CTRL"),
                     CellType = names(cellCountsInMS))

(Xsq2 <- chisq.test(t(CD8_subset_cellFrequenyData)))  # Prints test summary
Xsq2$observed   
Xsq2$expected  
Xsq2$residuals 
Xsq2$stdres 

stdResidualSig_cutoff <- function(nrow = 2, ncol = 2, alpha = .01){
  value <- qnorm(p=1-((alpha/2)/(nrow*ncol)))
  #output <- paste("Residual cutoff is", round(value, 3))
  return(round(value, 3))
}


pdf(file="cd8_df.js_chisquare_standardizedResiduals.pdf")
par(mar = c(4, 4.1, 4.1, 2.1))

barplot(sort(Xsq2$residuals[,"MS"]), las = 2,  ylim=c(-8,8),
        ylab="Standardized residual (MS)")

abline(h=c(stdResidualSig_cutoff(2,15),-stdResidualSig_cutoff(2,15)),col="red",lty=2)
axis(side=2, at=c(stdResidualSig_cutoff(2,15),-stdResidualSig_cutoff(2,15)), 
     labels=c(round(stdResidualSig_cutoff(2,15),2),round(-stdResidualSig_cutoff(2,15),2)),las=2,cex.axis=0.5)

dev.off()



# combined plot of cell type proportions per sample group
cellPropInMS = 100 * (sort(table(cd8_df.js@meta.data[cd8_df.js@meta.data$sampleGroup=="MS",]$seurat_clusters),decreasing=T)/sum(table(cd8_df.js@meta.data[cd8_df.js@meta.data$sampleGroup=="MS",]$seurat_clusters)))
cellPropInCTRL = 100 * (sort(table(cd8_df.js@meta.data[cd8_df.js@meta.data$sampleGroup=="CTRL",]$seurat_clusters),decreasing=T)/sum(table(cd8_df.js@meta.data[cd8_df.js@meta.data$sampleGroup=="CTRL",]$seurat_clusters)))

cellPropInMS_DF = data.frame(cellPropInMS)
colnames(cellPropInMS_DF) <- c("CellType","Relative_Frequency")
cellPropInMS_DF$sampleGroup="MS"

cellPropInCTRL_DF = data.frame(cellPropInCTRL)
colnames(cellPropInCTRL_DF) <- c("CellType","Relative_Frequency")
cellPropInCTRL_DF$sampleGroup="CTRL"

cellTypePropCombinedDF = rbind(cellPropInMS_DF,cellPropInCTRL_DF)



ggplot(data=cellTypePropCombinedDF, aes(x=CellType, y=Relative_Frequency, fill=sampleGroup)) +
  geom_bar(stat="identity", position=position_dodge()) + 
  #geom_point(aes(color= sampleGroup, shape=sampleGroup), alpha = 0.9, size=1.5,position = position_jitterdodge(jitter.width = 0.1)) + 
  scale_fill_manual(values=c("lightblue","Brown")) + 
  theme(legend.position = "top") +  theme_minimal() + theme(axis.text.x = element_text(angle = 90, hjust=1)) +
  ylab("Relative Frequency of CD8 subclusters (%)") + xlab("")

ggsave("cd8_df.js_cellTypePropCombinedDF.pdf", height = 210, width = 320, units = "mm")


#### Check the predicted TCR specificity profile in each cluster of CD8s ####


cd8_df.js@meta.data$predictedSpecificityVDJdbP <- NA
cd8_df.js@meta.data$predictedEpitopeGeneVDJdbP <- NA

cd8_df.js@meta.data$predictedSpecificityVDJdbP = msFNA.merged.tcr@meta.data$predictedSpecificityVDJdbP[match(rownames(cd8_df.js@meta.data),rownames(msFNA.merged.tcr@meta.data))]
cd8_df.js@meta.data$predictedEpitopeGeneVDJdbP = msFNA.merged.tcr@meta.data$predictedEpitopeGeneVDJdbP[match(rownames(cd8_df.js@meta.data),rownames(msFNA.merged.tcr@meta.data))]

cd8_df.js@meta.data$CTstrict = msFNA.merged.tcr@meta.data$CTstrict[match(rownames(cd8_df.js@meta.data),rownames(msFNA.merged.tcr@meta.data))]



TCRspecificityPredictions_byCluster = table(cd8_df.js@meta.data$predictedSpecificityVDJdbP,cd8_df.js@meta.data$seurat_clusters)

TCRspecificityPredictions_byClusterP = apply(TCRspecificityPredictions_byCluster,2,function(x) x/sum(x))

pheatmap::pheatmap(scale(TCRspecificityPredictions_byCluster),scale="row")




nHitsPerPathTypeList_cd8_df.js <- list()

for(cType in unique(cd8_df.js$seurat_clusters)){
  
  cd8_df.js.Temp <- cd8_df.js@meta.data[cd8_df.js$seurat_clusters==cType,]
    
 
  if(nrow(cd8_df.js.Temp) == 0)
    next
  
  nHitsPerPathTypeSamples <- NULL
  # nums to compare For all hit types
  for(sm in unique(cd8_df.js.Temp$sampleName)){
    
    samResForcTypeSample <- cd8_df.js.Temp[cd8_df.js.Temp$sampleName==sm,]
    
    samResForcTypeSample$predictedSpecificityVDJdbP[is.na(samResForcTypeSample$predictedSpecificityVDJdbP)] <- "NA"
    
    numsampTCRs <- length(rownames(samResForcTypeSample))
    
    nTCRsPerSamp <- c()
    
    for(pat in unique(VDJdblatest$antigen.species)){
      
      samResForcType <- samResForcTypeSample[samResForcTypeSample$predictedSpecificityVDJdbP==pat,] 
      
      numberOfhits <- ifelse(nrow(samResForcType) > 0, length(rownames(samResForcType)), 0)
      
      nTCRsPerSamp <- c(nTCRsPerSamp,(numberOfhits/numsampTCRs) * 1000)
      
    }
    
    nHitsPerPathTypeSamples <- cbind(nHitsPerPathTypeSamples,nTCRsPerSamp)
    
  }
  
  
  
  ### add to list  
  rownames(nHitsPerPathTypeSamples) <- unique(VDJdblatest$antigen.species)
  colnames(nHitsPerPathTypeSamples) <- unique(cd8_df.js.Temp$sampleName)
  
  msSams = which(grepl("MS",colnames(nHitsPerPathTypeSamples)))
  ctrSams = which(grepl("CTRL",colnames(nHitsPerPathTypeSamples)))
  
  if(length(msSams) < 3 | length(ctrSams) < 3) # at least three samples per group for which there is hit
    next
  
  
  nHitsPerPathTypeSamples_Comp = t(apply(nHitsPerPathTypeSamples + 1,1,function(x) c(wilcox.test(x[msSams],x[ctrSams])$p.value,
                                                                                     log2(mean(x[msSams])/mean(x[ctrSams])),log2(median(x[msSams])/median(x[ctrSams])))))
  
  # add t test p values for those that do not have constant values in all samples (which have NAN p values in wilcox test)
  idxForTtest = which(apply(nHitsPerPathTypeSamples,1,var) > 0)
  if(length(idxForTtest)==0)
    next
  
  nHitsPerPathTypeSamplesT.test = t(apply(nHitsPerPathTypeSamples[idxForTtest,,drop=F] + 1,1,function(x) c(t.test(x[msSams],x[ctrSams])$p.value,
                                                                                                           log2(mean(x[msSams])/mean(x[ctrSams])),log2(median(x[msSams])/median(x[ctrSams])))))
  
  nHitsPerPathTypeSamples_Comp[idxForTtest,1] <- nHitsPerPathTypeSamplesT.test[,1]
  
  
  
  colnames(nHitsPerPathTypeSamples_Comp) <- c("unpaired two-samples t.test.p.value","log2FC(meanMS/meanCTRL)","log2FC(medianMS/medianCTRL)")
  
  nHitsPerPathTypeSamplesFull = cbind(nHitsPerPathTypeSamples,nHitsPerPathTypeSamples_Comp)
  
  
  nHitsPerPathTypeList_cd8_df.js[[cType]] <- nHitsPerPathTypeSamplesFull
  
  
  
}

nHitsPerPathTypeList_cd8_df.js[["6"]] # CD8 cluster that shows the highest increase in frequency (in MS) shows significant increase in EBV specificity, and decrease to CMV specificity
nHitsPerPathTypeList_cd8_df.js[["5"]] # cluster with the second highest frequency increase in MS shows no differential specificity
nHitsPerPathTypeList_cd8_df.js[["9"]] # cluster with the highest exhaustion score shows no differential TCR specificity between MS and CTRls

nHitsPerPathTypeList_cd8_df.js[["11"]] # cluster that decreases the most in MS patients (in frequency) shows no differential specificity

nHitsPerPathTypeList_cd8_df.js[["14"]] # cluster that decreases the most in MS patients (in frequency) shows no differential specificity



## similar analysis for specificity predictions from Jani using TCRgp

TCRgp_predictions = read.csv("CD8_predictions_edited.csv")

cd8_df.js@meta.data$predictedSpecificityTCRgp <- NA
cd8_df.js@meta.data$predictedEpitopeGeneTCRgp <- NA

cd8_df.js@meta.data$predictedEpitopeGeneTCRgp_detailed <- NA
cd8_df.js@meta.data$target <- NA


# find the predictions for each of the CD8 cells based on CTstring and sample name
predIndice = sapply(1:nrow(cd8_df.js@meta.data),function(cl){
  
  cldata = cd8_df.js@meta.data[cl,]
  predIdx = which(TCRgp_predictions$CTstrict==cldata$CTstrict & TCRgp_predictions$sampleName==cldata$sampleName)
  
  if(length(predIdx) > 0){
    cd8_df.js@meta.data[cl,]$predictedSpecificityTCRgp <<- TCRgp_predictions[predIdx,]$virus
    cd8_df.js@meta.data[cl,]$predictedEpitopeGeneTCRgp <<- TCRgp_predictions[predIdx,]$virus_type
    cd8_df.js@meta.data[cl,]$predictedEpitopeGeneTCRgp_detailed <<- paste0(TCRgp_predictions[predIdx,]$cdr3b_pred_epitopes,",",TCRgp_predictions[predIdx,]$cdr3ab_pred_epitopes)
    cd8_df.js@meta.data[cl,]$target <<- TCRgp_predictions[predIdx,]$target
    
  }
  
  predIdx
})
predIndice = unlist(predIndice)

# heatmap of overall TCRgp prediction between clusters
TCRspecificityPredictions_tcrgp_byCluster = table(cd8_df.js@meta.data$predictedSpecificityTCRgp,cd8_df.js@meta.data$seurat_clusters)

TCRspecificityPredictions_tcrgp_byClusterP = apply(TCRspecificityPredictions_tcrgp_byCluster,2,function(x) x/sum(x))

pdf("cd8_df.js_TCRspecificityPredictions_tcrgp_byClusterP.pdf")
pheatmap::pheatmap(TCRspecificityPredictions_tcrgp_byClusterP,scale="row")

dev.off()

pdf("cd8_df.js_TCRspecificityPredictions_byClusterP.pdf")
pheatmap::pheatmap(TCRspecificityPredictions_byClusterP,scale="row")
dev.off()




nHitsPerPathTypeList_tcrgp_cd8_df.js <- list()

evaluatedTargets = unique(TCRgp_predictions$virus)

for(cType in unique(cd8_df.js$seurat_clusters)){
  
  cd8_df.js.Temp <- cd8_df.js@meta.data[cd8_df.js$seurat_clusters==cType,]
  
  
  if(nrow(cd8_df.js.Temp) == 0)
    next
  
  nHitsPerPathTypeSamples <- NULL
  # nums to compare For all hit types
  for(sm in unique(cd8_df.js.Temp$sampleName)){
    
    samResForcTypeSample <- cd8_df.js.Temp[cd8_df.js.Temp$sampleName==sm,]
    
    samResForcTypeSample$predictedSpecificityVDJdbP[is.na(samResForcTypeSample$predictedSpecificityTCRgp)] <- "NA"
    
    numsampTCRs <- length(rownames(samResForcTypeSample))
    
    nTCRsPerSamp <- c()
    
    for(pat in evaluatedTargets){
      
      samResForcType <- samResForcTypeSample[samResForcTypeSample$predictedSpecificityTCRgp==pat,] 
      
      numberOfhits <- ifelse(nrow(samResForcType) > 0, length(rownames(samResForcType)), 0)
      
      nTCRsPerSamp <- c(nTCRsPerSamp,(numberOfhits/numsampTCRs) * 1000)
      
    }
    
    nHitsPerPathTypeSamples <- cbind(nHitsPerPathTypeSamples,nTCRsPerSamp)
    
  }
  
  
  
  ### add to list  
  rownames(nHitsPerPathTypeSamples) <- evaluatedTargets
  colnames(nHitsPerPathTypeSamples) <- unique(cd8_df.js.Temp$sampleName)
  
  msSams = which(grepl("MS",colnames(nHitsPerPathTypeSamples)))
  ctrSams = which(grepl("CTRL",colnames(nHitsPerPathTypeSamples)))
  
  if(length(msSams) < 3 | length(ctrSams) < 3) # at least three samples per group for which there is hit
    next
  
  
  nHitsPerPathTypeSamples_Comp = t(apply(nHitsPerPathTypeSamples + 1,1,function(x) c(wilcox.test(x[msSams],x[ctrSams])$p.value,
                                                                                     log2(mean(x[msSams])/mean(x[ctrSams])),log2(median(x[msSams])/median(x[ctrSams])))))
  
  # add t test p values for those that do not have constant values in all samples (which have NAN p values in wilcox test)
  idxForTtest = which(apply(nHitsPerPathTypeSamples,1,var) > 0)
  if(length(idxForTtest)==0)
    next
  
  nHitsPerPathTypeSamplesT.test = t(apply(nHitsPerPathTypeSamples[idxForTtest,,drop=F] + 1,1,function(x) c(t.test(x[msSams],x[ctrSams])$p.value,
                                                                                                           log2(mean(x[msSams])/mean(x[ctrSams])),log2(median(x[msSams])/median(x[ctrSams])))))
  
  nHitsPerPathTypeSamples_Comp[idxForTtest,1] <- nHitsPerPathTypeSamplesT.test[,1]
  
  
  
  colnames(nHitsPerPathTypeSamples_Comp) <- c("unpaired two-samples t.test.p.value","log2FC(meanMS/meanCTRL)","log2FC(medianMS/medianCTRL)")
  
  nHitsPerPathTypeSamplesFull = cbind(nHitsPerPathTypeSamples,nHitsPerPathTypeSamples_Comp)
  
  
  nHitsPerPathTypeList_tcrgp_cd8_df.js[[cType]] <- nHitsPerPathTypeSamplesFull
  
  
  
}

nHitsPerPathTypeList_tcrgp_cd8_df.js[["6"]] # CD8 cluster that shows the highest increase in frequency (in MS) shows significant increase in EBV specificity, and decrease to CMV specificity
nHitsPerPathTypeList_tcrgp_cd8_df.js[["5"]] # cluster with the second highest frequency increase in MS shows no differential specificity
nHitsPerPathTypeList_tcrgp_cd8_df.js[["9"]] # cluster with the highest exhaustion score shows no differential TCR specificity between MS and CTRls

nHitsPerPathTypeList_tcrgp_cd8_df.js[["11"]] # cluster that decreases the most in MS patients (in frequency) shows no differential specificity



nHitsPerPathTypeList_cd8_df.js[["6"]]
nHitsPerPathTypeList_tcrgp_cd8_df.js[["6"]] 


nHitsPerPathTypeList_cd8_df.js[["0"]]
nHitsPerPathTypeList_tcrgp_cd8_df.js[["0"]] 

nHitsPerPathTypeList_tcrgp_cd8_df.js[["14"]] 


# check the distribution of cells of MS011 predicted to be EBV specific across the subclusters of CD8
# That is, what percentage of the EBV predicted cells in the sample belong to each CD8 sub cluster ? 

ebvTargetingCells_byCD8clustes = NULL
ebvTargetingCells_byCD8clustes_tcrgp = NULL

for(smp in unique(cd8_df.js$sampleName)){
    
  cd8_df.js.smp = subset(cd8_df.js, subset = sampleName == smp)
  ebvpredicted_insmp_cd8cls = table(cd8_df.js.smp$predictedSpecificityVDJdbP,cd8_df.js.smp$seurat_clusters)
  tcrgp_ebvpredicted_insmp_cd8cls = table(cd8_df.js.smp$predictedSpecificityTCRgp,cd8_df.js.smp$seurat_clusters)
  
  # what percentage of the EBV predicted cells in the sample belong to each CD8 sub cluster
  
  ebvpredicted_insmp_cd8cls_p = 100 * (ebvpredicted_insmp_cd8cls["EBV",] / sum(ebvpredicted_insmp_cd8cls["EBV",]))
  tcrgp_ebvpredicted_insmp_cd8cls_p = 100 * (tcrgp_ebvpredicted_insmp_cd8cls["EBV",] / sum(tcrgp_ebvpredicted_insmp_cd8cls["EBV",]))
  
  ebvTargetingCells_byCD8clustes = rbind(ebvTargetingCells_byCD8clustes,ebvpredicted_insmp_cd8cls_p)
  ebvTargetingCells_byCD8clustes_tcrgp = rbind(ebvTargetingCells_byCD8clustes_tcrgp,tcrgp_ebvpredicted_insmp_cd8cls_p)
  
  
}


rownames(ebvTargetingCells_byCD8clustes) <- unique(cd8_df.js$sampleName)
rownames(ebvTargetingCells_byCD8clustes_tcrgp) <- unique(cd8_df.js$sampleName)

ebvTargetingCells_byCD8clustes <- ebvTargetingCells_byCD8clustes[order(rownames(ebvTargetingCells_byCD8clustes)),]
ebvTargetingCells_byCD8clustes_tcrgp <- ebvTargetingCells_byCD8clustes_tcrgp[order(rownames(ebvTargetingCells_byCD8clustes_tcrgp)),]


pdf(paste0("cd8_df.js_cd8subclusters_allocationOfEbvTargettingCells_tcrdist.pdf"))
pheatmap::pheatmap(ebvTargetingCells_byCD8clustes, scale="none",cluster_rows = T, cluster_cols = F)
dev.off()

pdf(paste0("cd8_df.js_cd8subclusters_allocationOfEbvTargettingCells_tcrgp.pdf"))
pheatmap::pheatmap(ebvTargetingCells_byCD8clustes_tcrgp, scale="none",cluster_rows = T, cluster_cols = F)
dev.off()



## CD8 subcluster DE between MS and CTRLs ####

resdir = "cd8/"

table(cd8_df.js@meta.data$sampleGroup,cd8_df.js@meta.data$seurat_clusters)

for(cellT in unique(cd8_df.js@meta.data$seurat_clusters)){
  
  print(cellT)
  cd8_df.js_ct <- subset(x = cd8_df.js, subset = seurat_clusters == cellT)
  DefaultAssay(cd8_df.js_ct) <- "RNA"
  
  # single cell DE
  cd8_df.js_ct_DE = FindMarkers(cd8_df.js_ct,ident.1="MS",ident.2="CTRL",group.by = "sampleGroup")
  cd8_df.js_ct_DE_sig = cd8_df.js_ct_DE[cd8_df.js_ct_DE$p_val_adj < 0.05,]
  cd8_df.js_ct_DE_sig <- cd8_df.js_ct_DE_sig[order(cd8_df.js_ct_DE_sig$avg_log2FC,decreasing=T),]
  
  print(head(cd8_df.js_ct_DE_sig))
  
  write.csv(cd8_df.js_ct_DE,file=paste0(resdir,"CD8_js_subcluster_",cellT,"_MSvsCTRL_singleCellDE.csv"))
  write.csv(cd8_df.js_ct_DE_sig,file=paste0(resdir,"CD8_js_subcluster_",cellT,"_MSvsCTRL_singleCellDE_sig.csv"))
  
} 




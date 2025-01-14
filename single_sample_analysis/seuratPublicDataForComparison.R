
#### For comparison, we see the distributions of detected genes, UMIs, etc. in the example PBMC data used in seurat (https://satijalab.org/seurat/articles/pbmc3k_tutorial.html).

#pbmc.data <- Read10X(data.dir = "/scratch/project_2002480/SeuratAnalysis/publicPBMC/")
pbmc.data <- readRDS(file="/scratch/project_2005497/Dawit_scripts_MSFNA/seurat.pbmc.data")
numGenesDetectedPerCell.pbmc <- apply(pbmc.data,2, function(x) sum(x>0))

nfeaturesRaw_percentile <- ecdf(numGenesDetectedPerCell.pbmc)
nfeatures_min = 200 
nfeaturesRaw_max = nfeaturesRaw_percentile(200)
nfeaturesRaw_max = nfeaturesRaw_percentile(2500)



pbmc <- CreateSeuratObject(counts = pbmc.data, project = "pbmc3k", min.cells = 3, min.features = 200)
  
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")

mt_percentile <- ecdf(pbmc@meta.data$percent.mt)
mt_maxQuantile = round(mt_percentile(5),digits=2) # the typical max mt gene expression allowed in the pbmc data

# we do the same and estimate the min and max quantiles for nFeatures (genes). For now filtering by genes is enough, and filtering by UMI count not used.
nfeatures_percentile <- ecdf(pbmc@meta.data$nFeature_RNA)
#hist(pbmc@meta.data$nFeature_RNA)
nfeatures_min = nfeatures_percentile(200) 
nfeatures_max = nfeatures_percentile(2500)





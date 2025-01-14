library(Seurat)
library(ggplot2)
library(plyr)
library(ggpubr)
library(dplyr)


library(scRepertoire)


# load the merged Seurat object of the entire data ####

msFNAmerged <- readRDS(file = "data/msFNAmerged.rds")


# read the filtered BCR files into a list for BCR analysis #### 

# Use the directory where the sc BCR files are found for all samples (example naming: MS002_filtered_bcontig_annotations.csv)

filteredContigFileNames = list.files("/scRNAseq/vdjdata/",pattern ="_bcontig_",full.names=T,recursive = TRUE)



MS5_BCRcontigList <- lapply(filteredContigFileNames, function(x) read.csv(file = x,header=T))
names(MS5_BCRcontigList) <- sapply(filteredContigFileNames,
                                   function(y) paste0(unlist(strsplit(basename(y),"_"))[1],collapse="_"))


# # add sample barcode prefixes to the BCR data 

for(x in 1:length(MS5_BCRcontigList)){
  MS5_BCRcontigList[[x]]$barcode <- paste(names(MS5_BCRcontigList[x]),MS5_BCRcontigList[[x]]$barcode,sep="_")
}


# Now data is ready for analysis with scRepertoire:

# Combining the Contigs (of alpha and beta chain sequences for each cell are combined)

sampleGrp <- msFNAmerged@meta.data$sampleGroup[match(names(MS5_BCRcontigList),msFNAmerged@meta.data$orig.ident)]
sampleNms <- msFNAmerged@meta.data$sampleName[match(names(MS5_BCRcontigList),msFNAmerged@meta.data$orig.ident)]


combined_MS5_BCRcontigList <- combineBCR(MS5_BCRcontigList, 
                                         samples = sampleNms, 
                                         ID = sampleGrp, removeNA=T,removeMulti=T)


# modify the barcodes to match the barcodes in the seurat object here.
for(x in 1:length(combined_MS5_BCRcontigList)){
  combined_MS5_BCRcontigList[[x]]$barcode <- sapply(combined_MS5_BCRcontigList[[x]]$barcode,
                                                    function(x) paste0(unlist(strsplit(x,"_"))[c(-2,-3)],collapse="_"))
}
names(combined_MS5_BCRcontigList) <- sampleNms


## Combine BCR data to expression data
# Samples "MS001","CTRL012","CTRL013","CTRL014" are removed from analyses because they contained significantly less number of cells compared to other samples

samnamesSelected = samnames[!samnames %in% c("MS001","CTRL012","CTRL013","CTRL014")]
msFNAmerged.forBcr <- subset(x = msFNAmerged, subset = sampleName %in% samnamesSelected)

msFNA.merged.bcr <- combineExpression(combined_MS5_BCRcontigList, 
                                      msFNAmerged.forBcr, 
                                      cloneCall="CTstrict", 
                                      groupBy = "sample")

# Overall detection level of BCRs

# We have BCRs for 11048 of the 11462 immMain predicted B cells (96%)
perOfBellsWithBCRs = sum(!is.na(msFNA.merged.bcr@meta.data$barcode)) / table(msFNA.merged.bcr@meta.data$immMain)["B cells"]


cellsWithBCR <- msFNA.merged.bcr@meta.data$barcode[!is.na(msFNA.merged.bcr@meta.data$barcode)]

# entire B cells
msFNA.merged.bcrData = subset(msFNA.merged.bcr,subset = barcode %in% cellsWithBCR)

# Naive B cells
# 4495 of 4765 naives have BCR (94%)
msFNA.merged.bcrData.naive = subset(msFNA.merged.bcrData,subset = rnaclustersvsAdtManAnn2MaxAnn2 == "Naive B cells")

# Memory B cells
# 5382 of 6287 have BCR (lower BCR detection, 85%), our Mem B cell cluster is likely not pure, as expected
msFNA.merged.bcrData.mem = subset(msFNA.merged.bcrData,subset = rnaclustersvsAdtManAnn2MaxAnn2 == "Mem B cells")

# 662 of 686 plasmablasts have bcr (96%)
msFNA.merged.bcrData.pl = subset(msFNA.merged.bcrData,subset = rnaclustersvsAdtManAnn2MaxAnn2 == "Plasmablasts")

# 378 of 611 GC b cells have bcr (this is lower, 61%)
msFNA.merged.bcrData.gc = subset(msFNA.merged.bcrData,subset = rnaclustersvsAdtManAnn2MaxAnn2 == "GC B cells")



#### Save important data objects ####
saveRDS(combined_MS5_BCRcontigList, file = "data/combined_MS5_BCRcontigList.rds")
saveRDS(msFNA.merged.bcr, file = "data/msFNA.merged.bcr.rds")


### Complete BCR analysis for all B cell cell types ####


cellTypesWithBCRs <- names(which(table(msFNA.merged.bcrData@meta.data$rnaclustersvsAdtManAnn2MaxAnn2) > 100))
ResultDir = "BCRanalysis_Results_plotsWithSampleNames/"

normalized_shannon_index=function(f){
  
  tot_sum=sum(f)
  fi=f/tot_sum
  
  first <- fi * log(fi)
  sh <- abs(sum(first))
  #sh_effective=exp(sh) # effective size of repertoire
  nsh <- sh/log(length(f)) # normalizing entropy to number of unique TCRs

  return(nsh)
}

for(cbcr in c("allBCR",cellTypesWithBCRs)){
  
  print(cbcr)
  if(cbcr == "allBCR"){
    msFNA.merged.bcrTempData = subset(msFNA.merged.bcrData,subset = rnaclustersvsAdtManAnn2MaxAnn2 %in% cellTypesWithBCRs)
  }else{
    msFNA.merged.bcrTempData = subset(msFNA.merged.bcrData,subset = rnaclustersvsAdtManAnn2MaxAnn2 == cbcr)
  }
  combined_MS5_BCRcontigList_bcrTempData <- expression2List(msFNA.merged.bcrTempData,group="sampleName")
  
  write.to = paste0(ResultDir,paste0(cbcr,"_finalAnnotation"),"/")
  dir.create(write.to)
  
  # % of unique clonotypes 
  quantContig(combined_MS5_BCRcontigList_bcrTempData, cloneCall="CTstrict", scale = T)
  ggsave(paste0(write.to,cbcr,"_percentOfUniqueClonotypes.pdf"), width = 2, height = 4)
  
  # % of unique clonotypes by group
  quantContig(combined_MS5_BCRcontigList_bcrTempData, cloneCall="CTstrict", scale = T,group="sampleGroup")
  ggsave(paste0(write.to,cbcr,"_percentOfUniqueClonotypes_ByGroup.pdf"),  width = 2, height = 4)
  
  percOfuniqd = quantContig(combined_MS5_BCRcontigList_bcrTempData, cloneCall="CTstrict", scale = T,exportTable = T)
  percOfuniqd.ttest = t.test(percOfuniqd$scaled[grepl("MS",percOfuniqd$values)],percOfuniqd$scaled[grepl("CTRL",percOfuniqd$values)])
  
  write.csv(percOfuniqd,paste0(write.to,cbcr,"_percentOfUniqueClonotypes.csv"))
  write.csv(percOfuniqd.ttest$p.value,paste0(write.to,cbcr,"_percentOfUniqueClonotypes_ttestPval.csv"))
  
  percOfuniqd$sampleGroup <- ifelse(grepl("MS",percOfuniqd$values),"MS","CTRL")
  percOfuniqd$sampleGroup <- as.factor(percOfuniqd$sampleGroup)
  
  
  percOfuniqd_statp  = tryCatch({compare_means(scaled ~ sampleGroup, percOfuniqd,method="t.test",p.adjust.method="BH")},
                                error=function(e) {
                                  return(NA)
                                })
  
  if(class(percOfuniqd_statp)[1] == "tbl_df"){
    stat.test <- percOfuniqd_statp %>% mutate(y.position = rep(max(percOfuniqd$scaled) + 0.5,1))
    
    
    g <- ggplot(data = percOfuniqd, aes(x = sampleGroup,y = scaled,fill=sampleGroup))  
    g +  geom_boxplot(outlier.colour = NA) + 
      geom_text(label=percOfuniqd$values,size = 1,colour = "red",check_overlap = TRUE, position=position_jitter(width=0.50)) +
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      #scale_fill_manual(values=c("lightblue","Brown")) +
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      # geom_point(aes(shape = factor(percOfuniqd$values)),position=position_jitterdodge(jitter.width=0)) +
      theme_classic() + theme(legend.position = "none") + 
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Percent of Unique BCR Clonotypes (AA)") + ggtitle(cbcr) + NoLegend() + 
      
      # significance indicated with stars
      #stat_pvalue_manual(stat.test, x = "timePoint", hide.ns = TRUE,label = "p.signif",label.size=6,inherit.aes = FALSE,tip.length = 0)
      
      # significance with pvalues
      stat_pvalue_manual(stat.test, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)
    
    ggsave(paste0(write.to,cbcr,"_PercentageOfUniqueAAClonotypes_ForPaper.pdf"), width = 2, height = 4)
    
    
  }else{
    
    
    g <- ggplot(data = percOfuniqd, aes(x = sampleGroup,y = scaled,fill=sampleGroup))  
    g +  geom_boxplot(outlier.colour = NA) + geom_text(label=percOfuniqd$values,size = 1,colour = "red",check_overlap = TRUE, position=position_jitter(width=0.50)) +
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      #scale_fill_manual(values=c("lightblue","Brown")) +
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") +
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Percent of Unique BCR Clonotypes (AA)") + ggtitle(cbcr) 
    
    #ggsave(paste0(write.to,cbcr,"_PercentageOfUniqueAAClonotypes_ForPaper.pdf"), height = 210, width = 320, units = "mm")
    ggsave(paste0(write.to,cbcr,"_PercentageOfUniqueAAClonotypes_ForPaper.pdf"), width = 2, height = 4)
  }
  
  
  # number of GC BCRs species out of total BCR species (clonotypes) in each sample
  
  # get normalized shannon index values
  totalBCRperSample = sapply(percOfuniqd$values,function(x){
    
    msFNA.merged.bcrData.x = msFNA.merged.bcrData@meta.data[msFNA.merged.bcrData@meta.data$sampleName==x,]
    totalbcr= length(table(msFNA.merged.bcrData.x$CTstrict))
    totalbcr
  })
  
  percOfuniqd$totatBCRFromTotalBCRperc =  100 * (percOfuniqd$contigs/totalBCRperSample)
  
  percOfuniqd3.ttest = t.test(percOfuniqd$totatBCRFromTotalBCRperc[grepl("MS",percOfuniqd$values)],percOfuniqd$totatBCRFromTotalBCRperc[grepl("CTRL",percOfuniqd$values)])
  
  write.csv(percOfuniqd,paste0(write.to,cbcr,"_percentOfUniqueClonotypesAndRichness.csv"))
  
  
  # plot the BCR richness
  
  BCRrichness_statp  = tryCatch({compare_means(totatBCRFromTotalBCRperc ~ sampleGroup, percOfuniqd,method="t.test",p.adjust.method="BH")},
                                       error=function(e) {
                                         return(NA)
                                       })
  
  if(class(BCRrichness_statp)[1] == "tbl_df"){
    
    stat.test <- BCRrichness_statp %>% mutate(y.position = rep(max(percOfuniqd$totatBCRFromTotalBCRperc) + 0.001,1))
    
    g <- ggplot(data = percOfuniqd, aes(x = sampleGroup,y = totatBCRFromTotalBCRperc,fill=sampleGroup))  
    g +  geom_boxplot(outlier.colour = NA) + geom_text(label=percOfuniqd$values,size = 1,colour = "red",check_overlap = TRUE, position=position_jitter(width=0.50)) +
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      #scale_fill_manual(values=c("lightblue","Brown")) +
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      # geom_label(
      #   label=percOfuniqd$values,label.padding = unit(0.08, "lines"), 
      #   nudge_x = 0.30,nudge_y = -0.05,
      #   check_overlap = T,size=1
      # ) +
      theme_classic() + theme(legend.position = "none") +
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Richness (in percent of total BCR clonotypes in sample)") + ggtitle(cbcr) + 
      
      # significance indicated with stars
      #stat_pvalue_manual(stat.test, x = "timePoint", hide.ns = TRUE,label = "p.signif",label.size=6,inherit.aes = FALSE,tip.length = 0)
      
      # significance with pvalues
      stat_pvalue_manual(stat.test, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)
    #ggsave(paste0(write.to,cbcr,"_bcrdiversity_richness_forPaper.pdf"), height = 210, width = 320, units = "mm")
    
    ggsave(paste0(write.to,cbcr,"_bcrdiversity_richness_forPaper.pdf"), width = 2, height = 4)
    
    
  }else{
    
    g <- ggplot(data = percOfuniqd, aes(x = sampleGroup,y = totatBCRFromTotalBCRperc,fill=sampleGroup))  
    g +  geom_boxplot(outlier.colour = NA) + geom_text(label=percOfuniqd$values,size = 1,colour = "red",check_overlap = TRUE, position=position_jitter(width=0.50)) +
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      #scale_fill_manual(values=c("lightblue","Brown")) +
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") +
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Richness (in percent of total BCR clonotypes in sample)") + ggtitle(cbcr)
    
    #ggsave(paste0(write.to,cbcr,"_bcrdiversity_richness_forPaper.pdf"), height = 210, width = 320, units = "mm")
    
    ggsave(paste0(write.to,cbcr,"_bcrdiversity_richness_forPaper.pdf"), width = 2, height = 4)
    
    
  }
  
  # evaluate diversity
  # The shannon diversity is not normalized, needs to be normalized perhaps
  
  clonalDiversity(combined_MS5_BCRcontigList_bcrTempData,cloneCall = "CTstrict",group="sampleGroup",exportTable = F)
  #ggsave(paste0(write.to,cbcr,"_BCRdiversity.pdf"), height = 210, width = 320, units = "mm")
  ggsave(paste0(write.to,cbcr,"_BCRdiversity.pdf"), width = 2, height = 4)
  
  diversityBCRClones = clonalDiversity(combined_MS5_BCRcontigList_bcrTempData,cloneCall = "CTstrict",group="sampleGroup",exportTable = T)
  diversityBCRClones$sample=rownames(diversityBCRClones)
  
  diversityBCRClones.ttest = t.test(diversityBCRClones$Shannon[grepl("MS",diversityBCRClones$sample)],diversityBCRClones$Shannon[grepl("CTRL",diversityBCRClones$sample)])
  
  # get normalized shannon index values
  normShannonIndex = sapply(combined_MS5_BCRcontigList_bcrTempData,function(x){
    ftoPass= table(x$CTstrict)
    normalized_shannon_index(ftoPass)
    })

  diversityBCRClones$nShannon = normShannonIndex
  
  diversityBCRClones2.ttest = t.test(diversityBCRClones$nShannon[grepl("MS",diversityBCRClones$sample)],diversityBCRClones$nShannon[grepl("CTRL",diversityBCRClones$sample)])
  
  
  ggplot(diversityBCRClones,aes(x=sampleGroup,y=Shannon,color=sampleGroup)) + 
    geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + theme_bw() + labs(y="Diversity (Shannon)")
  
  ggsave(paste0(write.to,cbcr,"_BCRdiversity_Shannon.pdf"), height = 210, width = 320, units = "mm")
  #ggsave(paste0(write.to,cbcr,"_BCRdiversity_Shannon.pdf"), width = 2, height = 4)
  
  #
  ggplot(diversityBCRClones,aes(x=sampleGroup,y=Inv.Simpson,color=sampleGroup)) + 
    geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + theme_bw() + labs(y="Diversity (Inv.Simpson)")
  
  #ggsave(paste0(write.to,cbcr,"_BCRdiversity_Inv.Simpson.pdf"), height = 210, width = 320, units = "mm")
  
  ggsave(paste0(write.to,cbcr,"_BCRdiversity_Inv.Simpson.pdf"), width = 2, height = 4)
  
  
  write.csv(diversityBCRClones,paste0(write.to,cbcr,"_BCRdiversity.csv"))
  write.csv(diversityBCRClones.ttest$p.value,paste0(write.to,cbcr,"__BCRdiversity_ttestPval.csv"))
  
  # plot the normalized diversity for paper
  
  diversityBCRClones_statp  = tryCatch({compare_means(nShannon ~ sampleGroup, diversityBCRClones,method="t.test",p.adjust.method="BH")},
                                       error=function(e) {
                                         return(NA)
                                       })
  
  if(class(diversityBCRClones_statp)[1] == "tbl_df"){
    
    stat.test <- diversityBCRClones_statp %>% mutate(y.position = rep(max(diversityBCRClones$nShannon) + 0.001,1))
    
    g <- ggplot(data = diversityBCRClones, aes(x = sampleGroup,y = nShannon,fill=sampleGroup))  
    g +  geom_boxplot(outlier.colour = NA) + geom_text(label=diversityBCRClones$sample,size = 1,colour = "red",check_overlap = TRUE, position=position_jitter(width=0.50)) +
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      #scale_fill_manual(values=c("lightblue","Brown")) +
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") +
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Normalized Shannon diversity") + ggtitle(cbcr) + 
      
      # significance indicated with stars
      #stat_pvalue_manual(stat.test, x = "timePoint", hide.ns = TRUE,label = "p.signif",label.size=6,inherit.aes = FALSE,tip.length = 0)
      
      # significance with pvalues
      stat_pvalue_manual(stat.test, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)
    #ggsave(paste0(write.to,cbcr,"_bcrdiversity_normalizedShannon_forPaper.pdf"), height = 210, width = 320, units = "mm")
    
    ggsave(paste0(write.to,cbcr,"_bcrdiversity_normalizedShannon_forPaper.pdf"),  width = 2, height = 4)
    
    
  }else{
    
    g <- ggplot(data = diversityBCRClones, aes(x = sampleGroup,y = nShannon,fill=sampleGroup))  
    g +  geom_boxplot(outlier.colour = NA) + geom_text(label=diversityBCRClones$sample,size = 1,colour = "red",check_overlap = TRUE, position=position_jitter(width=0.50)) +
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      #scale_fill_manual(values=c("lightblue","Brown")) +
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") +
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Normalized Shannon diversity") + ggtitle(cbcr)
    
    #ggsave(paste0(write.to,cbcr,"_bcrdiversity_normalizedShannon_forPaper.pdf"), height = 210, width = 320, units = "mm")
    
    ggsave(paste0(write.to,cbcr,"_bcrdiversity_normalizedShannon_forPaper.pdf"), width = 2, height = 4)
    
  }
  
  
  
   # what kind of clones are in the BCR repertoire of GC etc.
  # These are calculated within the subset, and do not use the frequencies calculated earlier from the total BCR repertoire
  # for example for GC, even though they are diverse, they are so few that the space appears to be occupied by high expanded clones
  # That is also the case for controls
  clonalHomeostasis(combined_MS5_BCRcontigList_bcrTempData, cloneCall = "CTstrict",exportTable = F)
  ggsave(paste0(write.to,cbcr,"_BCRClonalHomeostasis.pdf"),width = 2, height = 4)
  
  clonalHomeostasisBCRClones = as.data.frame(clonalHomeostasis(combined_MS5_BCRcontigList_bcrTempData, cloneCall = "CTstrict",exportTable = T))
  
  colnames(clonalHomeostasisBCRClones) <- c("Rare","Small","Medium","Large","Hyperexpanded")
  
  clonalHomeostasisBCRClones$sample = rownames(diversityBCRClones)
  clonalHomeostasisBCRClones$sampleGroup = sapply(rownames(clonalHomeostasisBCRClones),function(x) ifelse(grepl("MS",x),"MS","CTRL"))
  
  ggplot(clonalHomeostasisBCRClones,aes(x=sampleGroup,y=Hyperexpanded,color=sampleGroup)) + 
    geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + theme_bw() + labs(y="Hyperexpanded (Relative Abundance)")
  
  ggsave(paste0(write.to,cbcr,"_BCRclonalHomeostatisHyperexpanded_byGroup.pdf"), width = 2, height = 4)
  
  ggplot(clonalHomeostasisBCRClones,aes(x=sampleGroup,y=Large,color=sampleGroup)) + 
    geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + theme_bw() + labs(y="Large (Relative Abundance)")
  
  ggsave(paste0(write.to,cbcr,"_BCRclonalHomeostatisLarge_byGroup.pdf"),  width = 2, height = 4)
  
  ggplot(clonalHomeostasisBCRClones,aes(x=sampleGroup,y=Medium,color=sampleGroup)) + 
    geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + theme_bw() + labs(y="Medium (Relative Abundance)")
  
  ggsave(paste0(write.to,cbcr,"_BCRclonalHomeostatisMedium_byGroup.pdf"), width = 2, height = 4)
  
  
  clonalHomeostasisBCRClonesHyper_test = t.test(clonalHomeostasisBCRClones$Hyperexpanded[clonalHomeostasisBCRClones$sampleGroup=="MS"],
                                           clonalHomeostasisBCRClones$Hyperexpanded[clonalHomeostasisBCRClones$sampleGroup=="CTRL"])

  clonalHomeostasisBCRClonesLarge_test = t.test(clonalHomeostasisBCRClones$Large[clonalHomeostasisBCRClones$sampleGroup=="MS"],
                                           clonalHomeostasisBCRClones$Large[clonalHomeostasisBCRClones$sampleGroup=="CTRL"])
  
  write.csv(clonalHomeostasisBCRClones,paste0(write.to,cbcr,"_clonalHomeostasisBCRClones.csv"))
  write.csv(clonalHomeostasisBCRClonesHyper_test$p.value,paste0(write.to,cbcr,"_clonalHomeostasisBCRClonesHyper_test.csv"))
  write.csv(clonalHomeostasisBCRClonesLarge_test$p.value,paste0(write.to,cbcr,"_clonalHomeostasisBCRClonesLarge_test.csv"))
  
  
}



### Percent clonality comparison (similar definition to lanz et al. 2022) ####

# https://www.nature.com/articles/s41586-022-04432-7#Sec11

LanzLike_clEstimatationData = msFNA.merged.bcrData@meta.data

LanzLike_clEstimatationData$LClonotypeNumber = NA
LanzLike_clEstimatationData$ClonotypeFrequency = NA

# identify Lanz et al like clonotypes across all samples and assign them clonal names
# sharing the same HC and LC V and J genes and exhibiting >70% amino acid identity within the HC and LC CDR3s
ClonotypeCount=0

for(i in 1:nrow(LanzLike_clEstimatationData)){
  
  print(i)
  if(!is.na(LanzLike_clEstimatationData[i,]$LClonotypeNumber))
    next
  
  currentBCR = LanzLike_clEstimatationData[i,]
  
  ctgenes = unlist(strsplit(currentBCR$CTgene,"[\\._]"))
  ctaas = unlist(strsplit(currentBCR$CTaa,"_"))
  
  currentBCR_HCV = ctgenes[grepl("IGHV",ctgenes)]
  currentBCR_HCJ = ctgenes[grepl("IGHJ",ctgenes)]
  currentBCR_LCV = ctgenes[grepl("IGKV",ctgenes) | grepl("IGLV",ctgenes)]
  currentBCR_LCJ = ctgenes[grepl("IGKJ",ctgenes) | grepl("IGLJ",ctgenes)]
  HC_CDR3aa = ctaas[1]
  LC_CDR3aa = ctaas[2]
  

  # The rest to compare to
  
  if(sum(is.na(LanzLike_clEstimatationData$LClonotypeNumber)) == 0)
    break
  
  LanzLike_clEstimatationData_nonAssignedToClonotypes = LanzLike_clEstimatationData[is.na(LanzLike_clEstimatationData$LClonotypeNumber),]
  

  HCVs = c()
  HCJs = c()
  LCVs = c()
  LCJs = c()
  HC_CDR3aas = c()
  LC_CDR3aas = c()
  
  for(m in 1:nrow(LanzLike_clEstimatationData_nonAssignedToClonotypes)){
      x = LanzLike_clEstimatationData_nonAssignedToClonotypes[m,]
      genesv = unlist(strsplit(x$CTgene,"[\\._]"))
      
      HCVs = c(HCVs,genesv[grepl("IGHV",genesv)])
      HCJs = c(HCJs,genesv[grepl("IGHJ",genesv)])
      LCVs = c(LCVs,genesv[grepl("IGKV",ctgenes) | grepl("IGLV",ctgenes)])
      LCJs = c(LCJs,genesv[grepl("IGKJ",ctgenes) | grepl("IGLJ",ctgenes)])
      
      aas = unlist(strsplit(x$CTaa,"_"))
    
      HC_CDR3aas = c(HC_CDR3aas,aas[1])
      LC_CDR3aas = c(LC_CDR3aas,aas[2])

    }
  
    similar_HCV = currentBCR_HCV == HCVs
    similar_HCJ = currentBCR_HCJ == HCJs
    similar_LCV = currentBCR_LCV == LCVs
    similar_LCJ = currentBCR_LCJ == LCJs
    
    
    # normalized levenstien distance between the current BCR and others to compared to
    # normalized by the maximum possible distance between two strings (results LD between 0-1)
    # then to find similarity we subtract from 1, 1- LD
    HCCDR3aaSims = 1 - stringdist::stringdist(HC_CDR3aa,HC_CDR3aas,method='lv') / sapply(HC_CDR3aas, function(x) max(nchar(HC_CDR3aa),nchar(x)))
    LCCDR3aaSims = 1 - stringdist::stringdist(LC_CDR3aa,LC_CDR3aas,method='lv') / sapply(LC_CDR3aas, function(x) max(nchar(LC_CDR3aa),nchar(x)))
    
    similar_HC_CDR3aa = HCCDR3aaSims > 0.7 
    similar_LC_CDR3aa = LCCDR3aaSims > 0.7
    
    currentBCR_clonotype_idx = similar_HCV & similar_HCJ & similar_LCV & similar_LCJ & similar_HC_CDR3aa & similar_LC_CDR3aa
    
    if(sum(currentBCR_clonotype_idx) > 1){
      
      # this means there are other BCRs similar to the current BCR with our requirements. And thus we can define them all as a clonotype
      # if only 1 BCR satisfies the requirement, it is the current BCR itself, so it is a single BCR and does not belong to a clonotype
      clonotypeBarcodes = rownames(LanzLike_clEstimatationData_nonAssignedToClonotypes[currentBCR_clonotype_idx,])
      clonotypeSize = length(clonotypeBarcodes)
      ClonotypeCount = ClonotypeCount + 1
      clonotypeLabel = paste0("Clonotype",ClonotypeCount)
      
      LanzLike_clEstimatationData[clonotypeBarcodes,]$LClonotypeNumber <- clonotypeLabel
      
    }else{
      clonotypeBarcodes = rownames(LanzLike_clEstimatationData_nonAssignedToClonotypes[currentBCR_clonotype_idx,])
      
      LanzLike_clEstimatationData[clonotypeBarcodes,]$LClonotypeNumber <- "single"
    }
    
  
  
}

# CLonotypes are identified with the lanz method, now calculate percent clonality and compare between groups

PercentClonalityLikeLanz <- NULL

for(s in unique(LanzLike_clEstimatationData$sampleName)){
  
  # all B cells
  
  LanzLike_clEstimatationData_allBcells = LanzLike_clEstimatationData[LanzLike_clEstimatationData$rnaclustersvsAdtManAnn2MaxAnn2 %in% cellTypesWithBCRs,]
  
  LanzLike_clEstimatationDataTemp = LanzLike_clEstimatationData_allBcells[LanzLike_clEstimatationData_allBcells$sampleName==s,]
  samName = s
  samGroup= unique(LanzLike_clEstimatationDataTemp$sampleGroup)
  cellTypes = "B cells"
  perClonality = (sum(LanzLike_clEstimatationDataTemp$LClonotypeNumber != "single") / nrow(LanzLike_clEstimatationDataTemp)) * 100  
  
  PercentClonalityLikeLanz = rbind(PercentClonalityLikeLanz,
                                   c(samName,samGroup,cellTypes,perClonality))
  
  # "Naive B cells"
  
  LanzLike_clEstimatationData_allBcells = LanzLike_clEstimatationData[LanzLike_clEstimatationData$rnaclustersvsAdtManAnn2MaxAnn2 == "Naive B cells",]
  
  LanzLike_clEstimatationDataTemp = LanzLike_clEstimatationData_allBcells[LanzLike_clEstimatationData_allBcells$sampleName==s,]
  samName = s
  samGroup= unique(LanzLike_clEstimatationDataTemp$sampleGroup)
  cellTypes = "Naive B cells"
  perClonality = (sum(LanzLike_clEstimatationDataTemp$LClonotypeNumber != "single") / nrow(LanzLike_clEstimatationDataTemp)) * 100  
  
  PercentClonalityLikeLanz = rbind(PercentClonalityLikeLanz,
                                   c(samName,samGroup,cellTypes,perClonality))
  
  
  # "Mem B cells"
  
  LanzLike_clEstimatationData_allBcells = LanzLike_clEstimatationData[LanzLike_clEstimatationData$rnaclustersvsAdtManAnn2MaxAnn2 == "Mem B cells",]
  
  LanzLike_clEstimatationDataTemp = LanzLike_clEstimatationData_allBcells[LanzLike_clEstimatationData_allBcells$sampleName==s,]
  samName = s
  samGroup= unique(LanzLike_clEstimatationDataTemp$sampleGroup)
  cellTypes = "Mem B cells"
  perClonality = (sum(LanzLike_clEstimatationDataTemp$LClonotypeNumber != "single") / nrow(LanzLike_clEstimatationDataTemp)) * 100  
  
  PercentClonalityLikeLanz = rbind(PercentClonalityLikeLanz,
                                   c(samName,samGroup,cellTypes,perClonality))
  
  # "Plasmablasts"
  
  LanzLike_clEstimatationData_allBcells = LanzLike_clEstimatationData[LanzLike_clEstimatationData$rnaclustersvsAdtManAnn2MaxAnn2 == "Plasmablasts",]
  
  LanzLike_clEstimatationDataTemp = LanzLike_clEstimatationData_allBcells[LanzLike_clEstimatationData_allBcells$sampleName==s,]
  samName = s
  samGroup= unique(LanzLike_clEstimatationDataTemp$sampleGroup)
  cellTypes = "Plasmablasts"
  perClonality = (sum(LanzLike_clEstimatationDataTemp$LClonotypeNumber != "single") / nrow(LanzLike_clEstimatationDataTemp)) * 100  
  
  PercentClonalityLikeLanz = rbind(PercentClonalityLikeLanz,
                                   c(samName,samGroup,cellTypes,perClonality))
  
  # "GC B cells" 
  
  LanzLike_clEstimatationData_allBcells = LanzLike_clEstimatationData[LanzLike_clEstimatationData$rnaclustersvsAdtManAnn2MaxAnn2 == "GC B cells",]
  
  LanzLike_clEstimatationDataTemp = LanzLike_clEstimatationData_allBcells[LanzLike_clEstimatationData_allBcells$sampleName==s,]
  samName = s
  samGroup= unique(LanzLike_clEstimatationDataTemp$sampleGroup)
  cellTypes = "GC B cells"
  perClonality = (sum(LanzLike_clEstimatationDataTemp$LClonotypeNumber != "single") / nrow(LanzLike_clEstimatationDataTemp)) * 100  
  
  PercentClonalityLikeLanz = rbind(PercentClonalityLikeLanz,
                                   c(samName,samGroup,cellTypes,perClonality))
  
  
  
}

colnames(PercentClonalityLikeLanz) <- c("samName","samGroup","cellTypes","Clonality")

PercentClonalityLikeLanzDF <- as.data.frame(PercentClonalityLikeLanz)

# testing

# all Bs
PercentClonalityLikeLanzDFtemp <- PercentClonalityLikeLanzDF[PercentClonalityLikeLanzDF$cellTypes=="B cells",]
percClonanity_AllBs = t.test(as.numeric(PercentClonalityLikeLanzDFtemp$Clonality[PercentClonalityLikeLanzDFtemp$samGroup=="MS"]),
       as.numeric(PercentClonalityLikeLanzDFtemp$Clonality[PercentClonalityLikeLanzDFtemp$samGroup=="CTRL"]))


# Naive B cells
PercentClonalityLikeLanzDFtemp <- PercentClonalityLikeLanzDF[PercentClonalityLikeLanzDF$cellTypes=="Naive B cells",]
percClonanity_NaiveBs = t.test(as.numeric(PercentClonalityLikeLanzDFtemp$Clonality[PercentClonalityLikeLanzDFtemp$samGroup=="MS"]),
                             as.numeric(PercentClonalityLikeLanzDFtemp$Clonality[PercentClonalityLikeLanzDFtemp$samGroup=="CTRL"]))

# Mem B cells
PercentClonalityLikeLanzDFtemp <- PercentClonalityLikeLanzDF[PercentClonalityLikeLanzDF$cellTypes=="Mem B cells",]
percClonanity_MemBs = t.test(as.numeric(PercentClonalityLikeLanzDFtemp$Clonality[PercentClonalityLikeLanzDFtemp$samGroup=="MS"]),
                               as.numeric(PercentClonalityLikeLanzDFtemp$Clonality[PercentClonalityLikeLanzDFtemp$samGroup=="CTRL"]))



# Plasmablasts
PercentClonalityLikeLanzDFtemp <- PercentClonalityLikeLanzDF[PercentClonalityLikeLanzDF$cellTypes=="Plasmablasts",]
percClonanity_Plasmablasts = t.test(as.numeric(PercentClonalityLikeLanzDFtemp$Clonality[PercentClonalityLikeLanzDFtemp$samGroup=="MS"]),
                             as.numeric(PercentClonalityLikeLanzDFtemp$Clonality[PercentClonalityLikeLanzDFtemp$samGroup=="CTRL"]))


# GC B cells
PercentClonalityLikeLanzDFtemp <- PercentClonalityLikeLanzDF[PercentClonalityLikeLanzDF$cellTypes=="GC B cells",]
percClonanity_GCBcells = t.test(as.numeric(PercentClonalityLikeLanzDFtemp$Clonality[PercentClonalityLikeLanzDFtemp$samGroup=="MS"]),
                                    as.numeric(PercentClonalityLikeLanzDFtemp$Clonality[PercentClonalityLikeLanzDFtemp$samGroup=="CTRL"]))


# plots
PercentClonalityLikeLanzDF$Clonality <- as.numeric(PercentClonalityLikeLanzDF$Clonality)
PercentClonalityLikeLanzDF$samGroup <- as.factor(PercentClonalityLikeLanzDF$samGroup)
PercentClonalityLikeLanzDF$cellTypes <- as.factor(PercentClonalityLikeLanzDF$cellTypes)
PercentClonalityLikeLanzDF$pertOfUniques <- 100 - as.numeric(PercentClonalityLikeLanzDF$Clonality)
PercentClonalityLikeLanzDF$ms11OrNot = as.character(PercentClonalityLikeLanzDF$samName == "MS011")

write.csv(PercentClonalityLikeLanzDF[,!colnames(PercentClonalityLikeLanzDF) %in% ("ms11OrNot")],
          file = "BcellTypes_Lanz_clonality_finalAnnotation.csv")


for(cellType in unique(PercentClonalityLikeLanzDF$cellTypes)){
  
  PercentClonalityLikeLanzDFTemp = PercentClonalityLikeLanzDF[PercentClonalityLikeLanzDF$cellTypes==cellType,]
  
  lanzPercOfUniq_statp  = tryCatch({compare_means(pertOfUniques ~ samGroup, PercentClonalityLikeLanzDFTemp,method="t.test",p.adjust.method="BH")},
                                       error=function(e) {
                                         return(NA)
                                       })
  
  stat.test <- lanzPercOfUniq_statp %>% mutate(y.position = rep(max(PercentClonalityLikeLanzDFTemp$pertOfUniques) + 5,1))
  
  p <- ggplot(PercentClonalityLikeLanzDFTemp, aes(x=samGroup, y=pertOfUniques,fill=samGroup,label=samName)) + 
    geom_boxplot(outlier.colour = NA) +
    #scale_fill_manual(values=c("lightblue","Brown")) +
    geom_point(position=position_jitterdodge(jitter.width=0), pch=21,aes(shape=factor(ms11OrNot))) +
    scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
    
    geom_text(hjust=-0.1, vjust=0,size = 2) + 
  
    #scale_shape_manual(values=c("TRUE" = 17,"FALSE" = 21)) +
    
    #geom_dotplot(binaxis='y', stackdir = "center", position = "dodge",dotsize=0.5,alpha=0.4) +
    theme(legend.position = "top") +  theme_minimal() + 
    theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
    ylab("Percent of Unique Clonotypes (AA)") + theme_bw() + NoLegend() + 
  
    # significance with pvalues
    stat_pvalue_manual(stat.test, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)
  
 
  ggsave(paste0(cellType,"_Lanz_clonality_finalAnnotation_withSamnames.pdf"), width=2, height = 4)
  
  
}



### Gene usage comparisons ####

forGeneComparison = msFNA.merged.bcrData


forGeneComparison@meta.data$HCVs <- sapply(forGeneComparison@meta.data$CTgene,function(x){
  gs = unlist(strsplit(x,"[\\._]"))
  gs[grepl("IGHV",gs)]
})


forGeneComparison@meta.data$HCJs <- sapply(forGeneComparison@meta.data$CTgene,function(x){
  gs = unlist(strsplit(x,"[\\._]"))
  gs[grepl("IGHJ",gs)]
})




forGeneComparison@meta.data$LCVs <- sapply(forGeneComparison@meta.data$CTgene,function(x){
  gs = unlist(strsplit(x,"[\\._]"))
  gs[grepl("IGKV",gs) | grepl("IGLV",gs)]
})


forGeneComparison@meta.data$LCJs <- sapply(forGeneComparison@meta.data$CTgene,function(x){
  gs = unlist(strsplit(x,"[\\._]"))
  gs[grepl("IGKJ",gs) | grepl("IGLJ",gs)]
})


# We only compared IG heavy change V genes as was done on lanz et al, other genes like IGHJ, IGLV and J can be done.

forGeneComparison@meta.data$rnaclustersvsAdtManAnn2MaxAnn2[forGeneComparison@meta.data$rnaclustersvsAdtManAnn2MaxAnn2=="Mem B cells"] <- "Memory B cells"

# for heavy chains

for(cbcr in c("all",cellTypesWithBCRs)){
  
  if(cbcr=="all"){
    forGeneComparisonTempData = subset(forGeneComparison,subset = rnaclustersvsAdtManAnn2MaxAnn2 %in% cellTypesWithBCRs)
  }else{
    forGeneComparisonTempData = subset(forGeneComparison,subset = rnaclustersvsAdtManAnn2MaxAnn2 == cbcr)
  }
  
  #combined_MS5_BCRcontigList_bcrTempData <- expression2List(msFNA.merged.bcrTempData,group="sampleName")
  
  HVSsfreqs = forGeneComparisonTempData@meta.data %>%
    group_by(sampleName,HCVs) %>% dplyr::summarise(cnt = n()) %>% mutate(freq = 100 * (cnt / sum(cnt)))
  
  HVSsfreqs$sampleGroup = forGeneComparisonTempData@meta.data$sampleGroup[match(HVSsfreqs$sampleName,forGeneComparisonTempData@meta.data$sampleName)]
  
  geneUsageComparisonPvals = c()
  for(gn in unique(HVSsfreqs$HCVs)){
    gnFreqInCTRL = HVSsfreqs$freq[HVSsfreqs$HCVs==gn & HVSsfreqs$sampleGroup=="CTRL"]
    gnFreqInMS = HVSsfreqs$freq[HVSsfreqs$HCVs==gn & HVSsfreqs$sampleGroup=="MS"]
    
    if(length(gnFreqInCTRL) >=3 && length(gnFreqInMS) >=3 ){
      ttestRes = t.test(gnFreqInMS,gnFreqInCTRL)
      resp = ttestRes$p.value
      names(resp) <- gn
      geneUsageComparisonPvals <- c(geneUsageComparisonPvals,resp)
    }
  }
  
  sigGenes = geneUsageComparisonPvals[geneUsageComparisonPvals < 0.05]
  # only IGHV3-11 shows an an increased usage in patients and MS associated genes like IGHV4-59 are decreased in MS
  geneUsageComparisonPvalsAdj = p.adjust(geneUsageComparisonPvals, "BH")
  
  
  ggplot(HVSsfreqs, aes(x=HCVs, y=freq, color=sampleGroup)) + 
    geom_boxplot(aes(colour = sampleGroup),outlier.shape = NA) + 
    geom_point(aes(color= sampleGroup, shape=sampleGroup), alpha = 0.9, size=1.5,position = position_jitterdodge(jitter.width = 0.1)) + 
    scale_colour_manual(values=c("lightblue","Brown")) +
    theme(legend.position = "top") +  theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab("Relative Frequency")
  
  ggsave(paste0("IGHVgeneusage_ttested_",cbcr,"_finalAnnotation.pdf"), height = 210, width = 320, units = "mm")
  
  
  sigGenesFreqData = HVSsfreqs[HVSsfreqs$HCVs %in% names(sigGenes),]
  write.csv(sigGenesFreqData,file=paste0("IGHVgeneusage_ttest_sigGenes_",cbcr,"_finalAnnotation.csv"))
  
 
  
}

# for light chains

for(cbcr in c("all",cellTypesWithBCRs)){
  
  if(cbcr=="all"){
    forGeneComparisonTempData = subset(forGeneComparison,subset = rnaclustersvsAdtManAnn2MaxAnn2 %in% cellTypesWithBCRs)
  }else{
    forGeneComparisonTempData = subset(forGeneComparison,subset = rnaclustersvsAdtManAnn2MaxAnn2 == cbcr)
  }
  
  #combined_MS5_BCRcontigList_bcrTempData <- expression2List(msFNA.merged.bcrTempData,group="sampleName")
  
  HVSsfreqs = forGeneComparisonTempData@meta.data %>%
    group_by(sampleName,LCVs) %>% dplyr::summarise(cnt = n()) %>% mutate(freq = 100 * (cnt / sum(cnt)))
  
  HVSsfreqs$sampleGroup = forGeneComparisonTempData@meta.data$sampleGroup[match(HVSsfreqs$sampleName,forGeneComparisonTempData@meta.data$sampleName)]
  
  geneUsageComparisonPvals = c()
  for(gn in unique(HVSsfreqs$LCVs)){
    gnFreqInCTRL = HVSsfreqs$freq[HVSsfreqs$LCVs==gn & HVSsfreqs$sampleGroup=="CTRL"]
    gnFreqInMS = HVSsfreqs$freq[HVSsfreqs$LCVs==gn & HVSsfreqs$sampleGroup=="MS"]
    
    if(length(gnFreqInCTRL) >=3 && length(gnFreqInMS) >=3 ){
      ttestRes = t.test(gnFreqInMS,gnFreqInCTRL)
      resp = ttestRes$p.value
      names(resp) <- gn
      geneUsageComparisonPvals <- c(geneUsageComparisonPvals,resp)
    }
  }
  
  sigGenes = geneUsageComparisonPvals[geneUsageComparisonPvals < 0.05]
  # only IGHV3-11 shows an an increased usage in patients and MS associated genes like IGHV4-59 are decreased in MS
  geneUsageComparisonPvalsAdj = p.adjust(geneUsageComparisonPvals, "BH")
  
  
  # plot IGHVgenes
  #ggplot(HVSsfreqs, aes(x=HCVs, y=freq, fill=sampleGroup)) + 
  # geom_boxplot() + coord_flip()
  #ggsave("/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/BCRanalysis_Results/Vgeneusage_allBcells.pdf", height = 210, width = 320, units = "mm")
  
  
  ggplot(HVSsfreqs, aes(x=LCVs, y=freq, color=sampleGroup)) + 
    geom_boxplot(aes(colour = sampleGroup),outlier.shape = NA) + 
    geom_point(aes(color= sampleGroup, shape=sampleGroup), alpha = 0.9, size=1.5,position = position_jitterdodge(jitter.width = 0.1)) + 
    scale_colour_manual(values=c("lightblue","Brown")) +
    theme(legend.position = "top") +  theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab("Relative Frequency")
  
  ggsave(paste0("IGLVgeneusage_ttested_",cbcr,"_finalAnnotation.pdf"), height = 210, width = 320, units = "mm")
  
  
  sigGenesFreqData = HVSsfreqs[HVSsfreqs$LCVs %in% names(sigGenes),]
  write.csv(sigGenesFreqData,file=paste0("IGLVgeneusage_ttest_sigGenes_",cbcr,"_finalAnnotation.csv"))
  
 
}


# for Heavy and light chains gene combinations
# there is only 8 such combinations that exist in at least 3 samples for each group, and that can be compared by t.test
forGeneComparison@meta.data$HC_LCgenes = paste(forGeneComparison@meta.data$HCVs,forGeneComparison@meta.data$HCJs,
                                               forGeneComparison@meta.data$LCVs,forGeneComparison@meta.data$LCJs,
                                               sep="_")

HC_LCgenes_counts = table(forGeneComparison@meta.data$HC_LCgenes,forGeneComparison@meta.data$sampleGroup)

head(HC_LCgenes_counts[order(HC_LCgenes_counts[,2],decreasing=T),],10)

head(HC_LCgenes_counts[order(HC_LCgenes_counts[,1],decreasing=T),],10)


for(cbcr in c("all",cellTypesWithBCRs)){
  
  if(cbcr=="all"){
    forGeneComparisonTempData = subset(forGeneComparison,subset = rnaclustersvsAdtManAnn2MaxAnn2 %in% cellTypesWithBCRs)
  }else{
    forGeneComparisonTempData = subset(forGeneComparison,subset = rnaclustersvsAdtManAnn2MaxAnn2 == cbcr)
  }
  
  #combined_MS5_BCRcontigList_bcrTempData <- expression2List(msFNA.merged.bcrTempData,group="sampleName")
  
  HVSsfreqs = forGeneComparisonTempData@meta.data %>%
    group_by(sampleName,HC_LCgenes) %>% dplyr::summarise(cnt = n()) %>% mutate(freq = 100 * (cnt / sum(cnt)))
  
  HVSsfreqs$sampleGroup = forGeneComparisonTempData@meta.data$sampleGroup[match(HVSsfreqs$sampleName,forGeneComparisonTempData@meta.data$sampleName)]
  
  geneUsageComparisonPvals = c()
  for(gn in unique(HVSsfreqs$HC_LCgenes)){
    gnFreqInCTRL = HVSsfreqs$freq[HVSsfreqs$HC_LCgenes==gn & HVSsfreqs$sampleGroup=="CTRL"]
    gnFreqInMS = HVSsfreqs$freq[HVSsfreqs$HC_LCgenes==gn & HVSsfreqs$sampleGroup=="MS"]
    
    if(length(gnFreqInCTRL) >=3 && length(gnFreqInMS) >=3 ){
      ttestRes = t.test(gnFreqInMS,gnFreqInCTRL)
      resp = ttestRes$p.value
      names(resp) <- gn
      geneUsageComparisonPvals <- c(geneUsageComparisonPvals,resp)
    }
  }
  
  sigGenes = geneUsageComparisonPvals[geneUsageComparisonPvals < 0.05]
  geneUsageComparisonPvalsAdj = p.adjust(geneUsageComparisonPvals, "BH")
  
 
  
  ggplot(HVSsfreqs, aes(x=LCVs, y=freq, color=sampleGroup)) + 
    geom_boxplot(aes(colour = sampleGroup),outlier.shape = NA) + 
    geom_point(aes(color= sampleGroup, shape=sampleGroup), alpha = 0.9, size=1.5,position = position_jitterdodge(jitter.width = 0.1)) + 
    scale_colour_manual(values=c("lightblue","Brown")) +
    theme(legend.position = "top") +  theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    ylab("Relative Frequency")
  
  ggsave(paste0("IGLVgeneusage_ttested_",cbcr,"_finalAnnotation.pdf"), height = 210, width = 320, units = "mm")
  
  
  sigGenesFreqData = HVSsfreqs[HVSsfreqs$LCVs %in% names(sigGenes),]
  write.csv(sigGenesFreqData,file=paste0("IGLVgeneusage_ttest_sigGenes_",cbcr,"_finalAnnotation.csv"))
  
 
}




# for Heavy and light chains gene combinations, without caring about the J genes
# there is only 8 such combinations that exist in at least 3 samples for each group, and that can be compared by t.test
forGeneComparison@meta.data$HC_LC_Vgenes = paste(forGeneComparison@meta.data$HCVs,
                                               forGeneComparison@meta.data$LCVs,
                                               sep="_")

HC_LC_Vgenes_counts = table(forGeneComparison@meta.data$HC_LC_Vgenes,forGeneComparison@meta.data$sampleGroup)

head(HC_LC_Vgenes_counts[order(HC_LC_Vgenes_counts[,2],decreasing=T),],30)

head(HC_LC_Vgenes_counts[order(HC_LC_Vgenes_counts[,1],decreasing=T),],10)


for(cbcr in c("all",cellTypesWithBCRs)){
  
  if(cbcr=="all"){
    forGeneComparisonTempData = subset(forGeneComparison,subset = rnaclustersvsAdtManAnn2MaxAnn2 %in% cellTypesWithBCRs)
  }else{
    forGeneComparisonTempData = subset(forGeneComparison,subset = rnaclustersvsAdtManAnn2MaxAnn2 == cbcr)
  }
  
  #combined_MS5_BCRcontigList_bcrTempData <- expression2List(msFNA.merged.bcrTempData,group="sampleName")
  
  HVSsfreqs = forGeneComparisonTempData@meta.data %>%
    group_by(sampleName,HC_LC_Vgenes) %>% dplyr::summarise(cnt = n()) %>% mutate(freq = 100 * (cnt / sum(cnt)))
  
  HVSsfreqs$sampleGroup = forGeneComparisonTempData@meta.data$sampleGroup[match(HVSsfreqs$sampleName,forGeneComparisonTempData@meta.data$sampleName)]
  
  geneUsageComparisonPvals = c()
  for(gn in unique(HVSsfreqs$HC_LC_Vgenes)){
    gnFreqInCTRL = HVSsfreqs$freq[HVSsfreqs$HC_LC_Vgenes==gn & HVSsfreqs$sampleGroup=="CTRL"]
    gnFreqInMS = HVSsfreqs$freq[HVSsfreqs$HC_LC_Vgenes==gn & HVSsfreqs$sampleGroup=="MS"]
    
    if(length(gnFreqInCTRL) >=3 && length(gnFreqInMS) >=3 ){
      ttestRes = t.test(gnFreqInMS,gnFreqInCTRL)
      resp = ttestRes$p.value
      names(resp) <- gn
      geneUsageComparisonPvals <- c(geneUsageComparisonPvals,resp)
    }
  }
  
  sigGenes = geneUsageComparisonPvals[geneUsageComparisonPvals < 0.05]
  geneUsageComparisonPvalsAdj = p.adjust(geneUsageComparisonPvals, "BH")
  
  ggplot(HVSsfreqs[HVSsfreqs$HC_LC_Vgenes %in% names(sigGenes),], aes(x=HC_LC_Vgenes, y=freq, color=sampleGroup)) + 
    geom_boxplot(aes(colour = sampleGroup),outlier.shape = NA) + 
    geom_point(aes(color= sampleGroup, shape=sampleGroup), alpha = 0.9, size=1.5,position = position_jitterdodge(jitter.width = 0.1)) + 
    scale_colour_manual(values=c("lightblue","Brown")) +
    theme(legend.position = "top") +  theme_minimal() + 
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1,size = 5)) + 
    ylab("Relative Frequency") +  coord_flip()
  
  ggsave(paste0("HC_LC_Vgenes_usage_ttested_",cbcr,"_finalAnnotation.pdf"), height = 210, width = 320, units = "mm")
  
  
  sigGenesFreqData = HVSsfreqs[HVSsfreqs$HC_LC_Vgenes %in% names(sigGenes),]
  write.csv(sigGenesFreqData,file=paste0("HC_LC_Vgenes_usage_ttested_sigGenes_",cbcr,"_finalAnnotation.csv"))
  
  
}




#### Analysis after manuscript peer-review ####

##1. Evaluate Lanz et al's GLialCAM specific antibody MS39p2w174 in our BCRs ##

Lanz_MS39p2w174_likeBCR = msFNA.merged.bcrData@meta.data

# We also do the analysis at these levels

msFNAmerged.goodSamples.gcb = readRDS("data/msFNAmerged.goodSamples.bcells4")

msFNAmerged.goodSamples.memNaivB = readRDS("data/msFNAmerged.goodSamples.membcells")

bcellsTouse = c(colnames(msFNAmerged.goodSamples.memNaivB),colnames(msFNAmerged.goodSamples.gcb))

Lanz_MS39p2w174_likeBCR_onlyIn_memNaiveB_gcb = Lanz_MS39p2w174_likeBCR[rownames(Lanz_MS39p2w174_likeBCR) %in% bcellsTouse, ]


# find clones that are the same as or in the same clonotype as MS39p2w174 (from Lanz et al) 

Lanz_MS39p2w174_HCV = "IGHV3-7"
Lanz_MS39p2w174_HCJ = "IGHJ4"
Lanz_MS39p2w174_LCV = "IGKV2-30"
Lanz_MS39p2w174_LCJ = "IGKJ1"

# from http://aligncdr.labshare.cn/aligncdr/cdrs.php and https://www.nature.com/articles/s41586-022-04432-7/figures/11
# the CDR3s have conserved C at the start in both chains, also all our CDR3s of BCRs in our data have conserved C and W at HC and LC respectively
# but T at the second position of HC CDR3 aa is the second common after A, it is more frequent than the others: table(sapply(HC_CDR3aas, function(mll) unlist(strsplit(mll,""))[2]))

# for LC CDr3aas, C is again the only AA at the beginning, and F is the most frequent at the end which MS39p2w174 also has: https://www.nature.com/articles/s41586-022-04432-7/figures/11
# 

Lanz_MS39p2w174_HC_CDR3aa = "CTRDPPYFDNW" 
Lanz_MS39p2w174_LC_CDR3aa = "CMQGSHWPVTF"


Lanz_MS39p2w174_likeBCR$MS39p2w174Clone = NA

# find similar clonotypes to MS39p2w174 

HCVs = c()
HCJs = c()
LCVs = c()
LCJs = c()
HC_CDR3aas = c()
LC_CDR3aas = c()

for(m in 1:nrow(Lanz_MS39p2w174_likeBCR)){
  x = Lanz_MS39p2w174_likeBCR[m,]
  genesv = unlist(strsplit(x$CTgene,"[\\._]"))
  
  HCVs = c(HCVs,genesv[grepl("IGHV",genesv)])
  HCJs = c(HCJs,genesv[grepl("IGHJ",genesv)])
  LCVs = c(LCVs,genesv[grepl("IGKV",ctgenes) | grepl("IGLV",ctgenes)])
  LCJs = c(LCJs,genesv[grepl("IGKJ",ctgenes) | grepl("IGLJ",ctgenes)])
  
  aas = unlist(strsplit(x$CTaa,"_"))
  
  HC_CDR3aas = c(HC_CDR3aas,aas[1])
  LC_CDR3aas = c(LC_CDR3aas,aas[2])
  
}

similar_HCV = Lanz_MS39p2w174_HCV == HCVs
similar_HCJ = Lanz_MS39p2w174_HCJ == HCJs
similar_LCV = Lanz_MS39p2w174_LCV == LCVs
similar_LCJ = Lanz_MS39p2w174_LCJ == LCJs

#table(sapply(HC_CDR3aas, function(mll) unlist(strsplit(mll,""))[1]))

# sum(Lanz_MS39p2w174_HC_CDR3aa == HC_CDR3aas)  # no exact matching clonotypes to MS39p2w174 heavy chain CDR3
# sum(Lanz_MS39p2w174_LC_CDR3aa == LC_CDR3aas)  # no exact matching clonotypes to MS39p2w174 light chain CDR3


# normalized levenstien distance between the current BCR and others to compared to
# normalized by the maximum possible distance between two strings (results LD between 0-1)
# then to find similarity we subtract from 1, 1- LD
HCCDR3aaSims = 1 - stringdist::stringdist(Lanz_MS39p2w174_HC_CDR3aa,HC_CDR3aas,method='lv') / sapply(HC_CDR3aas, function(x) max(nchar(Lanz_MS39p2w174_HC_CDR3aa),nchar(x)))
LCCDR3aaSims = 1 - stringdist::stringdist(Lanz_MS39p2w174_LC_CDR3aa,LC_CDR3aas,method='lv') / sapply(LC_CDR3aas, function(x) max(nchar(Lanz_MS39p2w174_LC_CDR3aa),nchar(x)))

# similarity at 95 percentile and above
# since we detected very few clonotypes at 0.7 similarity with CDR3s of the glialcam and ebna1 specific clonotype MS39p2w174, we decided to compare
# all clonotypes in our data with similarity scores above 95 percentile to MS39p2w174, and compare their frequencies between patient groups

simCutoff95Perc_HC = quantile(HCCDR3aaSims,seq(0,1,0.01))[96]
simCutoff95Perc_LC = quantile(LCCDR3aaSims,seq(0,1,0.01))[96]

similar_HC_CDR3aa = HCCDR3aaSims > simCutoff95Perc_HC 
similar_LC_CDR3aa = LCCDR3aaSims > simCutoff95Perc_LC

# Comparison of MS39p2w174 like repertoire using strict similarity measure (same HCV, HCJ, LCV, LCJ genes, and high similarity in HC and LC AA CDR3s in the top 5 percentile of similarities from 11k clonotypes)
# This appears to show no difference bewteen groups, in fact at higher cutoffs no similar clonotypes are detected in all samples.

  similarClonotype_idx_strict = similar_HCV & similar_HCJ & similar_LCV & similar_LCJ & similar_HC_CDR3aa & similar_LC_CDR3aa
  
  #table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_strict,]$sampleGroup)/table(Lanz_MS39p2w174_likeBCR$sampleGroup)

  
# Comparison of MS39p2w174 like repertoire using only same HCV, LCV genes (V genes of heavy and light genes), and high similarity in HC and LC AA CDR3s in the top 5 percentile of similarities from 11k clonotypes)
# This identifies only 2 clonotypes both from MS patients (1 each from MS004 in USM and MS005), both memory B cells 
  
  similarClonotype_idx_vGeneAndCDR3Only = similar_HCV & similar_LCV & similar_HC_CDR3aa & similar_LC_CDR3aa
  
  detectedSimilarClonotypes = rownames(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_vGeneAndCDR3Only,]) 
  
    msFNAmerged.goodSamples.memNaivB@meta.data[detectedSimilarClonotypes,]
    msFNAmerged.goodSamples.gcb@meta.data[detectedSimilarClonotypes,]

    
# Comparison of MS39p2w174 like repertoire using only same HCV (only heavy chain V gene), and high similarity in HC and LC AA CDR3s in the top 5 percentile of similarities from 11k clonotypes)
# This identifies 21 clonotypes from both MS and Controls but with no difference in their frequency between the patient groups
  
    
  similarClonotype_idx_heavyChainOnly = similar_HCV & similar_HC_CDR3aa 
  
  table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_heavyChainOnly,]$sampleName)/table(Lanz_MS39p2w174_likeBCR$sampleName)
  
  table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_heavyChainOnly,]$sampleGroup)/table(Lanz_MS39p2w174_likeBCR$sampleGroup)
  

# Comparison of MS39p2w174 like repertoire using only CDR3 sequences, high similarity in HC and LC AA CDR3s in the top 5 percentile of similarities from 11k clonotypes)
# This identifies 20 clonotypes from both patient groups, with slightly higher frequency in controls (also the same above) but no overall statistical difference
  
  
  similarClonotype_idx_CD3sOnly = similar_HC_CDR3aa & similar_LC_CDR3aa
  
  table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_CD3sOnly,]$sampleGroup)/table(Lanz_MS39p2w174_likeBCR$sampleGroup)
  

# Conclusion: No clones that were exact matches of MS39p2w174 were found in all patients. When using strict similarity criteria of same V and J gene usage for both H and L chains
# and 0.7 and above similarity (or in top 5% percentile similarity) for CDR3s of both chains, or just CDR3aas, or just heavy chain similarity etc as seen above, there was no difference
# in the frequency of similar clones to MS39p2w174 between patient groups. For the case where CDr3aa similarity above 0.7 was used, there was no detection of any similar clone in all samples and cases

# For other suggested auto Antigenes in MS, we couldn't find BCR sequences for autoantiboies yet for anoctamin 2

  
  

# similar_HC_CDR3aa = HCCDR3aaSims > 0.7 
# similar_LC_CDR3aa = LCCDR3aaSims > 0.7
# 
# # Comparison of MS39p2w174 like repertoire using strict similarity measure (same HCV, HCJ, LCV, LCJ genes, and high similarity in HC and LC AA CDR3s in the top 5 percentile of similarities from 11k clonotypes)
# # This appears to show no difference bewteen groups, in fact at higher cutoffs no similar clonotypes are detected in all samples.
# 
# similarClonotype_idx_strict = similar_HCV & similar_HCJ & similar_LCV & similar_LCJ & similar_HC_CDR3aa & similar_LC_CDR3aa
# 
# #table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_strict,]$sampleGroup)/table(Lanz_MS39p2w174_likeBCR$sampleGroup)
# 
# 
# # Comparison of MS39p2w174 like repertoire using only same HCV, LCV genes (V genes of heavy and light genes), and high similarity in HC and LC AA CDR3s in the top 5 percentile of similarities from 11k clonotypes)
# # This identifies only 2 clonotypes both from MS patients (1 each from MS004 in USM and MS005), both memory B cells 
# 
# similarClonotype_idx_vGeneAndCDR3Only = similar_HCV & similar_LCV & similar_HC_CDR3aa & similar_LC_CDR3aa
# 
# detectedSimilarClonotypes = rownames(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_vGeneAndCDR3Only,]) 
# 
# msFNAmerged.goodSamples.memNaivB@meta.data[detectedSimilarClonotypes,]
# msFNAmerged.goodSamples.gcb@meta.data[detectedSimilarClonotypes,]
# 
# 
# # Comparison of MS39p2w174 like repertoire using only same HCV (only heavy chain V gene), and high similarity in HC and LC AA CDR3s in the top 5 percentile of similarities from 11k clonotypes)
# # This identifies 21 clonotypes from both MS and Controls but with no difference in their frequency between the patient groups
# 
# 
# similarClonotype_idx_heavyChainOnly = similar_HCV & similar_HC_CDR3aa 
# 
# table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_heavyChainOnly,]$sampleName)/table(Lanz_MS39p2w174_likeBCR$sampleName)
# 
# table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_heavyChainOnly,]$sampleGroup)/table(Lanz_MS39p2w174_likeBCR$sampleGroup)
# 
# 
# # Comparison of MS39p2w174 like repertoire using only CDR3 sequences, high similarity in HC and LC AA CDR3s in the top 5 percentile of similarities from 11k clonotypes)
# # This identifies 20 clonotypes from both patient groups, with slightly higher frequency in controls (also the same above) but no overall statistical difference
# 
# 
# similarClonotype_idx_CD3sOnly = similar_HC_CDR3aa & similar_LC_CDR3aa
# 
# table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_CD3sOnly,]$sampleGroup)/table(Lanz_MS39p2w174_likeBCR$sampleGroup)
# 


# Using a less strict similarity cutoff for the CDR3aas (similarity in the top 25 percentiles)

simCutoff95Perc_HC = quantile(HCCDR3aaSims,seq(0,1,0.01))[76]
simCutoff95Perc_LC = quantile(LCCDR3aaSims,seq(0,1,0.01))[76]

similar_HC_CDR3aa = HCCDR3aaSims > simCutoff95Perc_HC 
similar_LC_CDR3aa = LCCDR3aaSims > simCutoff95Perc_LC

# Comparison of MS39p2w174 like repertoire using strict similarity measure (same HCV, HCJ, LCV, LCJ genes, and high similarity in HC and LC AA CDR3s in the top 25 percentile of similarities from 11k clonotypes)
# This appears to show no difference bewteen groups, just one clonotype in one MS and one CTRL samples

similarClonotype_idx_strict = similar_HCV & similar_HCJ & similar_LCV & similar_LCJ & similar_HC_CDR3aa & similar_LC_CDR3aa

table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_strict,]$sampleGroup)/table(Lanz_MS39p2w174_likeBCR$sampleGroup)


# Comparison of MS39p2w174 like repertoire using only same HCV, LCV genes (V genes of heavy and light genes), and high similarity in HC and LC AA CDR3s in the top 25 percentile of similarities from 11k clonotypes)
# This identifies only 5 clonotypes, 4 in MS and 1 in CTRL, 0.04 vs 0.03 % of total repertoires MS vs CTRls 
# 2 of clonotypes in MS are USM and 1 SM types (MS004, MS004, MS002), and 1 USM in CTRL015

similarClonotype_idx_vGeneAndCDR3Only = similar_HCV & similar_LCV & similar_HC_CDR3aa & similar_LC_CDR3aa

table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_vGeneAndCDR3Only,]$sampleGroup)/table(Lanz_MS39p2w174_likeBCR$sampleGroup)

detectedSimilarClonotypes = rownames(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_vGeneAndCDR3Only,]) 

msFNAmerged.goodSamples.memNaivB@meta.data[detectedSimilarClonotypes,]
msFNAmerged.goodSamples.gcb@meta.data[detectedSimilarClonotypes,]


# Comparison of MS39p2w174 like repertoire using only same HCV (only heavy chain V gene), and high similarity in HC AA CDR3s in the top 25 percentile of similarities from 11k clonotypes)
# This identifies 98 clonotypes from both MS and Controls, 0.8 vs 0.9 % in ctrls vs MS, but when looking at each individual, there is no difference in frequency


similarClonotype_idx_heavyChainOnly = similar_HCV & similar_HC_CDR3aa 

table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_heavyChainOnly,]$sampleName)/table(Lanz_MS39p2w174_likeBCR$sampleName)

table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_heavyChainOnly,]$sampleGroup)/table(Lanz_MS39p2w174_likeBCR$sampleGroup)


# Comparison of MS39p2w174 like repertoire using only CDR3 sequences, high similarity in HC and LC AA CDR3s in the top 25 percentile of similarities from 11k clonotypes)
# This identifies 495 clonotypes from both patient groups, with slightly higher frequency in controls (also the same above) but no overall statistical difference


similarClonotype_idx_CD3sOnly = similar_HC_CDR3aa & similar_LC_CDR3aa

table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_CD3sOnly,]$sampleName)/table(Lanz_MS39p2w174_likeBCR$sampleName)

table(Lanz_MS39p2w174_likeBCR[similarClonotype_idx_CD3sOnly,]$sampleGroup)/table(Lanz_MS39p2w174_likeBCR$sampleGroup)

# Conclusion: there is no difference in frequency of clonotypes similar to MS39p2w174 even when using a looser definition for similarity of the CDr3 sequences. Only slighly increased frequencies
# in MS when using both Vgene and CDR3s of the chains. Generally, the results show no statistical difference in the frequency of clonotypes that are similar to MS39p2w174 between MS and CTRLs
  
# Save analysis workspaces
# save.image("MS5FNA_autoAntibodySearch.RData")
# load("MS5FNA_autoAntibodySearch.RData")



#### confirm the b cell identity of the GCb subcluster of GC b cells.
msFNAmerged.goodSamples.bcells4.GCb = subset(msFNAmerged.goodSamples.gcb,subset = GC_type=="GC.B")

table(msFNAmerged.goodSamples.bcells4.GCb@meta.data$immFine)

# There are 109 GC.B cells from all samples
# Majority are suggested to be b cells using both immFine or immMain annotation
# 10 are tregs with immFine, 1 Th1, 1 Tfh, 8 dendritic cells, 24 unannotated, the rest are different B cell types

# 31 of the 109 have BCR data
msFNA.merged.bcrData.GCb = msFNA.merged.bcrData[,colnames(msFNAmerged.goodSamples.bcells4.GCb)]

write.csv(msFNA.merged.bcrData.GCb@meta.data,
          "GC.B_subsetcellsWithBCRs.csv")


#
msFNA.tcr.tcrs = readRDS("data/msFNA.merged.tcr.rds")

msFNA.tcr.tcrs.GCb = msFNA.tcr.tcrs[,colnames(msFNAmerged.goodSamples.bcells4.GCb)]

# 7 of the 109 GC.B cells had TCRs detected. For all the 7 except one  only the b chain tcr was detected and the alpha chain was not detected.
nrow(msFNA.tcr.tcrs.GCb@meta.data[!is.na(msFNA.tcr.tcrs.GCb@meta.data$barcode),])

write.csv(msFNA.tcr.tcrs.GCb@meta.data[!is.na(msFNA.tcr.tcrs.GCb@meta.data$barcode),],
          "GC.B_subsetcellsWithTCRs.csv")




##2. Check if the transcriptome of EBV-infected lymphoblastoid cells lines (LCLs) from https://elifesciences.org/articles/62586 is enriched in MS ##

# use B-cell subclustering data
msFNAmerged.goodSamples.gcb = readRDS("data/msFNAmerged.goodSamples.bcells4")

msFNAmerged.goodSamples.memNaivB = readRDS("data/msFNAmerged.goodSamples.membcells")

# then process each public LCL data using the script provided by the authors of paper https://elifesciences.org/articles/62586
# (/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/Bcell_comparisons_toPublicData/Source_Code_File_1.R)
# From here we get the upregulated markers of the lytic cluster: host_lyticClusterMarkers and 
# the host lytic genes after removing all EBV genes: lyticClusterMarkers_sortByFC_withoutEBVgenes

host_lytic_genes_main = c("SFN", "MIER2", "NFATC1", "NHLH1", "SGK1")

VlnPlot(msFNAmerged.goodSamples.gcb, features = host_lytic_genes_main, ncol = 4, pt.size=0.5)

VlnPlot(msFNAmerged.goodSamples.memNaivB, features = host_lytic_genes_main[5], pt.size=0.5,split.by = "sampleGroup")

# first check which clusters in our B cell subclustered dataset are enriched with host lytic markers
msFNAmerged.goodSamples.gcb.markers <- FindAllMarkers(msFNAmerged.goodSamples.gcb,only.pos = T,logfc.threshold = 0.1)
msFNAmerged.goodSamples.gcb.markers.sig <- msFNAmerged.goodSamples.gcb.markers[msFNAmerged.goodSamples.gcb.markers$p_val_adj < 0.05,]

msFNAmerged.goodSamples.gcb.markers.sig[msFNAmerged.goodSamples.gcb.markers.sig$gene %in% host_lytic_genes_main ,]

totalNumberOfGenes_inLCLsample = nrow(lcl777_b958.ccr)
numberOflyticGenesAsMarkersInLCLSample = length(host_lyticClusterMarkers) 

msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment = NULL

for(cls in unique(msFNAmerged.goodSamples.gcb.markers.sig$cluster)){
  

  cls_marker_genes = msFNAmerged.goodSamples.gcb.markers.sig$gene[msFNAmerged.goodSamples.gcb.markers.sig$cluster==cls]
  nHostLyticCells_inCluster = length(intersect(host_lyticClusterMarkers,cls_marker_genes))
  
  percLytic = nHostLyticCells_inCluster/length(cls_marker_genes)
  
  print(cls)
  print(percLytic)
  print(length(cls_marker_genes))
 
  
  # Fishers exact test to check the proportion of lytic genes in the total set of genes we have in the reference and in our cluster
  
  dat <- data.frame(
    "nonLytic_markers" = c(totalNumberOfGenes_inLCLsample-numberOflyticGenesAsMarkersInLCLSample, length(cls_marker_genes)-nHostLyticCells_inCluster),
    "Lytic_markers" = c(numberOflyticGenesAsMarkersInLCLSample, nHostLyticCells_inCluster),
    row.names = c("LCLsample", "Our_bGC_data_cls"),
    stringsAsFactors = FALSE
  )

  ebvLyticmarkers_Enrichment_test <- fisher.test(dat)
  
  print(ebvLyticmarkers_Enrichment_test)
  
  cls_res = c(cls,length(cls_marker_genes),nHostLyticCells_inCluster,percLytic,ebvLyticmarkers_Enrichment_test$p.value,ebvLyticmarkers_Enrichment_test$estimate,paste0(ebvLyticmarkers_Enrichment_test$conf.int,collapse=" - "))
  
  msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment = rbind(msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment,cls_res)
  
  
}

colnames(msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment) <- c("gcb.cluster","nMarkers","nLyticMarkers","proportionOfLyticMarkers","FishersExact_pvalue","oddsRatio","OddsRatio_Conf.interval")

# so even though LZ gc BC cluster is the only cluster that has one of the main host lytic genes as marker, SGK1. Enrichment
# analysis of lytic genes in each cluster of gc b cell cluster shows that the GC.B cluster is overwhelmingly enriched (36%) by lytic markers
# and it was the only gc B cell cluster that was more abundant in MS compared to controls.


msFNAmerged.goodSamples.gcb.clusters.MSvsCTRL_lyticMarkerEnrichment = NULL
msFNAmerged.goodSamples.gcb.clusters.CTRLvsMS_lyticMarkerEnrichment = NULL

for(cls in unique(msFNAmerged.goodSamples.gcb.markers.sig$cluster)){
  
  msFNAmerged.goodSamples.gcb.cls = subset(msFNAmerged.goodSamples.gcb, subset = GC_type==cls)
  
  msFNAmerged.goodSamples.gcb.cls.markers <- FindMarkers(msFNAmerged.goodSamples.gcb.cls,only.pos = T,logfc.threshold = 0.1,ident.1 = "MS",ident.2 = "CTRL",group.by = "sampleGroup")
  msFNAmerged.goodSamples.gcb.cls.markers.sig <- msFNAmerged.goodSamples.gcb.cls.markers[msFNAmerged.goodSamples.gcb.cls.markers$p_val < 0.05,]
  
  msFNAmerged.goodSamples.gcb.cls.markers.sig[rownames(msFNAmerged.goodSamples.gcb.cls.markers.sig) %in% host_lyticClusterMarkers,]
  
  
  cls_marker_genes = rownames(msFNAmerged.goodSamples.gcb.cls.markers.sig)
  nHostLyticCells_inCluster = length(intersect(host_lyticClusterMarkers,cls_marker_genes))
  
  percLytic = nHostLyticCells_inCluster/length(cls_marker_genes)
  
  print(cls)
  print(percLytic)
  print(length(cls_marker_genes))
  
  
  # Fishers exact test to check the proportion of lytic genes in the total set of genes we have in the reference and in our cluster
  
  dat <- data.frame(
    "nonLytic_markers" = c(totalNumberOfGenes_inLCLsample-numberOflyticGenesAsMarkersInLCLSample, length(cls_marker_genes)-nHostLyticCells_inCluster),
    "Lytic_markers" = c(numberOflyticGenesAsMarkersInLCLSample, nHostLyticCells_inCluster),
    row.names = c("LCLsample", "Our_bGC_data_cls"),
    stringsAsFactors = FALSE
  )
  
  ebvLyticmarkers_Enrichment_test <- fisher.test(dat)
  
  print(ebvLyticmarkers_Enrichment_test)
  
  cls_res = c(cls,length(cls_marker_genes),nHostLyticCells_inCluster,percLytic,ebvLyticmarkers_Enrichment_test$p.value,ebvLyticmarkers_Enrichment_test$estimate,paste0(ebvLyticmarkers_Enrichment_test$conf.int,collapse=" - "))
  
  msFNAmerged.goodSamples.gcb.clusters.MSvsCTRL_lyticMarkerEnrichment = rbind(msFNAmerged.goodSamples.gcb.clusters.MSvsCTRL_lyticMarkerEnrichment,cls_res)
  
  
}

colnames(msFNAmerged.goodSamples.gcb.clusters.MSvsCTRL_lyticMarkerEnrichment) <- c("gcb.cluster","nMarkers","nLyticMarkers","proportionOfLyticMarkers","FishersExact_pvalue","oddsRatio","OddsRatio_Conf.interval")




# plot to show the lytic markers in our clusters

writeDir = "Bcell_comparisons_toPublicData/"

msFNAmerged.goodSamples.gcbExp = AggregateExpression(msFNAmerged.goodSamples.gcb,assays = "RNA",slot="counts")

msFNAmerged.goodSamples.gcbExpNormalized = apply(msFNAmerged.goodSamples.gcbExp$RNA,2,function(x) 10000 * (x/sum(x)))

msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers <- msFNAmerged.goodSamples.gcbExpNormalized[rownames(msFNAmerged.goodSamples.gcbExpNormalized) %in% host_lyticClusterMarkers,]

msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers = msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers[rowMeans(msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers) > 0,]


msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers_scaled = t(scale(t(as.matrix(msFNAmerged.goodSamples.gcbExpNormalized))))


pdf(file = paste0(writeDir,"msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers_scaled_heatmap_aggregate.pdf"),width=6, height=10)

toplyticGenes = rownames(lyticClusterMarkers_sortByFC_withoutEBVgenes)[1:50]

ComplexHeatmap::pheatmap(msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers_scaled[rownames(msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers_scaled) %in% toplyticGenes,],
                         scale="none",fontsize_row = 5,fontsize_col=6,
                         cluster_rows = T,column_names_side = c("top"),
                         color =colorRampPalette(c("blue", "white", "red"))(100))

dev.off()


msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers2 = msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers[,!grepl("prolif",colnames(msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers))]

msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers2 = msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers2[rowMeans(msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers2) > 0,]


msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers2_scaled = t(scale(t(as.matrix(msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers2) + 1)))

msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers2_scaled_nonRibosomalgenes = msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers2_scaled[!grepl("RP",rownames(msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers2_scaled)),]

pdf(file = paste0(writeDir,"msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers2_scaled_heatmap_aggregate.pdf"),width=6, height=10)


ComplexHeatmap::pheatmap(msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers2_scaled,scale="none",fontsize_row = 6,fontsize_col=6,
                         cluster_rows = T,cluster_cols = T,column_names_side = c("bottom"),
                         color =colorRampPalette(c("blue", "white", "red"))(100))

dev.off()



# I decided to use the above plot (after excluding prolif.T cells)
# We show lytic all host lytic genes in our data and their expression in normalized sum aggregage expression by sample/cellsubset
# The lytic genes with highest FC appear to be highly expressed in MS subsets. And generally most lytic genes have higher expression in GC.B subcluster (of both MS and controls)


# but it might be good to score the clusters using the top lytic genes

toplyticGenes = lyticClusterMarkers_sortByFC_withoutEBVgenes$gene[lyticClusterMarkers_sortByFC_withoutEBVgenes$avg_log2FC >= 1]

Lytic_features <- list(HostLyticgenesTop=toplyticGenes,HostLyticgenesTop=host_lyticClusterMarkers)

toplyticGenes[toplyticGenes %in% rownames(msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers2_scaled)]

msFNAmerged.goodSamples.gcb <-  AddModuleScore(object = msFNAmerged.goodSamples.gcb,assay = "RNA",
                                                          features = Lytic_features,ctrl = 100,
                                                          name = 'lyticScore')



nLyticsPergcbCluster_gcb <- msFNAmerged.goodSamples.gcb@meta.data %>% group_by(sampleGroup,GC_type) %>% 
   dplyr::summarise(n=n(),nLytics = sum(lyticScore1 > 0), mainlyticscoreMean = mean(lyticScore1),alllyticscoreMean = mean(lyticScore2)) %>% mutate(pmainLytics = (nLytics/n) * 100) %>% arrange(desc(alllyticscoreMean))

write.csv(lyticClusterMarkers_sortByFC_withoutEBVgenes,
          file = paste0(writeDir,"LCL_lyticCell_hostlyticMarkers.csv"))

msFNAmerged.goodSamples.gcb[["sampleGroup_GCtype"]] <- paste(msFNAmerged.goodSamples.gcb@meta.data$sampleGroup,msFNAmerged.goodSamples.gcb@meta.data$GC_type,sep="_")

RidgePlot(subset(msFNAmerged.goodSamples.gcb,subset=GC_type != "prolif.T"), features = "lyticScore1",group.by=c("sampleGroup_GCtype")) + ggtitle("cell scoring of main host lytic markers")

ggsave(file = paste0(writeDir,"msFNAmerged.goodSamples.gcb_mainLyticmarkers_cellScore.pdf")) 


RidgePlot(subset(msFNAmerged.goodSamples.gcb,subset=GC_type != "prolif.T"), features = "lyticScore2") + ggtitle("cell scoring of all host lytic markers")
ggsave(file = paste0(writeDir,"msFNAmerged.goodSamples.gcb_allLyticmarkers_cellScore.pdf"))


RidgePlot(subset(msFNAmerged.goodSamples.gcb,subset=GC_type != "prolif.T"), features = "lyticScore2",group.by=c("sampleGroup_GCtype")) + ggtitle("cell scoring of all host lytic markers")
ggsave(file = paste0(writeDir,"msFNAmerged.goodSamples.gcb_allLyticmarkers_cellScore_sampleGroup_GCtype.pdf"))

# main lytic genes

msFNAmerged.goodSamples.gcb_clsAvgExpD_mainlyticMarkers <- msFNAmerged.goodSamples.gcbExpNormalized[rownames(msFNAmerged.goodSamples.gcbExpNormalized) %in% host_lytic_genes_main,]


msFNAmerged.goodSamples.gcb_clsAvgExpD_mainlyticMarkers_scaled = t(scale(t(as.matrix(msFNAmerged.goodSamples.gcb_clsAvgExpD_mainlyticMarkers) + 1)))



pdf(file = paste0(writeDir,"msFNAmerged.goodSamples.gcb_clsAvgExpD_mainlyticMarkers_scaled_heatmap_aggregate.pdf"),width=6, height=10)


ComplexHeatmap::pheatmap(msFNAmerged.goodSamples.gcb_clsAvgExpD_mainlyticMarkers_scaled,scale="none",fontsize_row = 6,fontsize_col=6,
                         cluster_rows = T,cluster_cols = F,column_names_side = c("bottom"),
                         color =colorRampPalette(c("blue", "white", "red"))(100))

dev.off()



# for memNaive B cells


msFNAmerged.goodSamples.memNaivBExp = AggregateExpression(msFNAmerged.goodSamples.memNaivB,group.by = c('sampleGroup', 'memnaive_type'),assays = "RNA",slot="counts")

msFNAmerged.goodSamples.memNaivBExpNormalized = apply(msFNAmerged.goodSamples.memNaivBExp$RNA,2,function(x) 10000 * (x/sum(x)))

msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers <- msFNAmerged.goodSamples.memNaivBExpNormalized[rownames(msFNAmerged.goodSamples.memNaivBExpNormalized) %in% host_lyticClusterMarkers,]

msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers = msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers[rowMeans(msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers) > 0,]


msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers_scaled = t(scale(t(as.matrix(msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers) + 1)))

#msFNAmerged.goodSamples.gcb_clsAvgExpD_lyticMarkers_scaled = t(scale(t(as.matrix(msFNAmerged.goodSamples.gcbExpNormalized))))


pdf(file = paste0(writeDir,"msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers_scaled_heatmap_aggregate.pdf"),width=6, height=10)

toplyticGenes = rownames(lyticClusterMarkers_sortByFC_withoutEBVgenes)[1:50]

ComplexHeatmap::pheatmap(msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers_scaled[rownames(msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers_scaled) %in% toplyticGenes,],
                         scale="none",fontsize_row = 5,fontsize_col=6,
                         cluster_rows = T,column_names_side = c("top"),
                         color =colorRampPalette(c("blue", "white", "red"))(100))

dev.off()


msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers2 = msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers[,!grepl("prolif",colnames(msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers))]

msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers2 = msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers2[rowMeans(msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers2) > 0,]


msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers2_scaled = t(scale(t(as.matrix(msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers2) + 1)))

msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers2_scaled_nonRibosomalgenes = msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers2_scaled[!grepl("RP",rownames(msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers2_scaled)),]

pdf(file = paste0(writeDir,"msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers2_scaled_heatmap_aggregate.pdf"),width=6, height=10)


ComplexHeatmap::pheatmap(msFNAmerged.goodSamples.memNaivB_clsAvgExpD_lyticMarkers2_scaled,scale="none",fontsize_row = 6,fontsize_col=6,
                         cluster_rows = T,cluster_cols = T,column_names_side = c("bottom"),
                         color =colorRampPalette(c("blue", "white", "red"))(100))

dev.off()



# lytic scoring for bcell memory naives

msFNAmerged.goodSamples.memNaivB <-  AddModuleScore(object = msFNAmerged.goodSamples.memNaivB,assay = "RNA",
                                               features = Lytic_features,ctrl = 100,
                                               name = 'lyticScore')

nLyticsPermemNaivBCluster_memNaivB <- msFNAmerged.goodSamples.memNaivB@meta.data %>% group_by(sampleGroup,memnaive_type) %>% 
  dplyr::summarise(n=n(),nLytics = sum(lyticScore1 > 0), mainlyticscoreMean = mean(lyticScore1),alllyticscoreMean = mean(lyticScore2)) %>% mutate(pmainLytics = (nLytics/n) * 100) %>% arrange(desc(alllyticscoreMean))


msFNAmerged.goodSamples.memNaivB[["sampleGroup_memnaivetype"]] <- paste(msFNAmerged.goodSamples.memNaivB@meta.data$sampleGroup,msFNAmerged.goodSamples.memNaivB@meta.data$memnaive_type,sep="_")

RidgePlot(msFNAmerged.goodSamples.memNaivB, features = "lyticScore1") + ggtitle("cell scoring of main host lytic markers")
ggsave(file = paste0(writeDir,"msFNAmerged.goodSamples.memNaivB_mainLyticmarkers_cellScore.pdf")) 


RidgePlot(msFNAmerged.goodSamples.memNaivB, features = "lyticScore2") + ggtitle("cell scoring of all host lytic markers")
ggsave(file = paste0(writeDir,"msFNAmerged.goodSamples.memNaivB_allLyticmarkers_cellScore.pdf"))


RidgePlot(msFNAmerged.goodSamples.memNaivB, features = "lyticScore2",group.by=c("sampleGroup_memnaivetype")) + ggtitle("cell scoring of all host lytic markers")
ggsave(file = paste0(writeDir,"msFNAmerged.goodSamples.memNaivB_allLyticmarkers_cellScore_sampleGroup_memnaivetype.pdf"))



## Host lytic gene enrichment analysis for all subclusters of naive and memory B cells: msFNAmerged.goodSamples.memNaivB

msFNAmerged.goodSamples.memNaivB.markers <- FindAllMarkers(msFNAmerged.goodSamples.memNaivB,only.pos = T,logfc.threshold = 0.1)
msFNAmerged.goodSamples.memNaivB.markers.sig <- msFNAmerged.goodSamples.memNaivB.markers[msFNAmerged.goodSamples.memNaivB.markers$p_val_adj < 0.05,]

msFNAmerged.goodSamples.memNaivB.markers.sig[msFNAmerged.goodSamples.memNaivB.markers.sig$gene %in% host_lytic_genes_main ,] # none of the main lytic genes were markers of the big B cell subsets



msFNAmerged.goodSamples.memNaivB.clusters.lyticMarkerEnrichment = NULL

for(cls in unique(msFNAmerged.goodSamples.memNaivB.markers.sig$cluster)){
  
  
  cls_marker_genes = msFNAmerged.goodSamples.memNaivB.markers.sig$gene[msFNAmerged.goodSamples.memNaivB.markers.sig$cluster==cls]
  nHostLyticCells_inCluster = length(intersect(host_lyticClusterMarkers,cls_marker_genes))
  
  percLytic = nHostLyticCells_inCluster/length(cls_marker_genes)
  
  print(cls)
  print(percLytic)
  print(length(cls_marker_genes))
  
  
  # Fishers exact test to check the proportion of lytic genes in the total set of genes we have in the reference and in our cluster
  
  dat <- data.frame(
    "nonLytic_markers" = c(totalNumberOfGenes_inLCLsample-numberOflyticGenesAsMarkersInLCLSample, length(cls_marker_genes)-nHostLyticCells_inCluster),
    "Lytic_markers" = c(numberOflyticGenesAsMarkersInLCLSample, nHostLyticCells_inCluster),
    row.names = c("LCLsample", "Our_bGC_data_cls"),
    stringsAsFactors = FALSE
  )
  
  ebvLyticmarkers_Enrichment_test <- fisher.test(dat)
  
  print(ebvLyticmarkers_Enrichment_test)
  
  cls_res = c(cls,length(cls_marker_genes),nHostLyticCells_inCluster,percLytic,ebvLyticmarkers_Enrichment_test$p.value,ebvLyticmarkers_Enrichment_test$estimate,paste0(ebvLyticmarkers_Enrichment_test$conf.int,collapse=" - "))
  
  msFNAmerged.goodSamples.memNaivB.clusters.lyticMarkerEnrichment = rbind(msFNAmerged.goodSamples.memNaivB.clusters.lyticMarkerEnrichment,cls_res)
  
  
}

colnames(msFNAmerged.goodSamples.memNaivB.clusters.lyticMarkerEnrichment) <- c("bcell.cluster","nMarkers","nLyticMarkers","proportionOfLyticMarkers","FishersExact_pvalue","oddsRatio","OddsRatio_Conf.interval")

## Just like the GC.B cluster GC cluster, the DN B cell cluster and to some extent the SM B clusters are characterized by over-representation of 
## host lytic markers -> suggesting that these are definitely more prevalent in MS and are in lytic stage. We don't know how specific the lytic genes are to EBV lytic stage or if they are shared with other
## virus lytic stages. But overall, host lytic genes are over-represented in b cell clusters that show significantly more abundance in MS patients compared to controls.

write.csv(msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment,
          "msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment.csv")


write.csv(msFNAmerged.goodSamples.memNaivB.clusters.lyticMarkerEnrichment,
          "msFNAmerged.goodSamples.memNaivB.clusters.lyticMarkerEnrichment.csv")




# plot host lytic gene enrichment

memNaivB.clusters.lyticMarkerEnrichment = as.data.frame(msFNAmerged.goodSamples.memNaivB.clusters.lyticMarkerEnrichment)
memNaivB.clusters.lyticMarkerEnrichment$Bcellstype = "Extrafollicular"
colnames(memNaivB.clusters.lyticMarkerEnrichment)[1] <- "cluster"

gcb.clusters.lyticMarkerEnrichment = as.data.frame(msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment)
gcb.clusters.lyticMarkerEnrichment$Bcellstype = "Follicular"
colnames(gcb.clusters.lyticMarkerEnrichment)[1] <- "cluster"

Bcell.clusters.lyticMarkerEnrichment = rbind(memNaivB.clusters.lyticMarkerEnrichment,gcb.clusters.lyticMarkerEnrichment)

Bcell.clusters.lyticMarkerEnrichment$oddsRatio = as.numeric(Bcell.clusters.lyticMarkerEnrichment$oddsRatio)
Bcell.clusters.lyticMarkerEnrichment$p.value = as.numeric(Bcell.clusters.lyticMarkerEnrichment$FishersExact_pvalue) < 0.05

ggplot(data = Bcell.clusters.lyticMarkerEnrichment[Bcell.clusters.lyticMarkerEnrichment$Bcellstype=="Extrafollicular",], aes(x = cluster, y = oddsRatio,fill=p.value)) + coord_flip() + 
  geom_point(shape = 21, colour = "black",size=4) + 
  geom_hline(yintercept=1, linetype="dashed", color = "grey") + 
  scale_fill_manual(values = c("white","red")) +
  theme_classic() + 
  ylab("OR") + guides(fill=guide_legend(title="p.value < 0.05")) +
  facet_grid(cols  = vars(Bcellstype))

ggsave("Extrafollicular_lyticmarker_overrepresentation.pdf")


ggplot(data = Bcell.clusters.lyticMarkerEnrichment[Bcell.clusters.lyticMarkerEnrichment$Bcellstype=="Follicular",], aes(x = cluster, y = oddsRatio,fill=p.value)) + coord_flip() + 
  geom_point(shape = 21, colour = "black",size=4) + 
  geom_hline(yintercept=1, linetype="dashed", color = "grey") + 
  scale_fill_manual(values = c("white","red")) +
  theme_classic() + 
  ylab("OR") + guides(fill=guide_legend(title="p.value < 0.05")) +
  facet_grid(cols  = vars(Bcellstype))

ggsave("Follicular_lyticmarker_overrepresentation.pdf")


write.csv(Bcell.clusters.lyticMarkerEnrichment,
          file = "Bcell.clusters.lyticMarkerEnrichment.csv")




## split the host lytic markers Odds ratio (OR) by sample or patient group.

# For this, since the markers of the clusters are defined for clusters (containing cells from multiple samples), marker detection
# is done within each sample separately. And enrichment OR calculated for each sample.


msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment.persample = NULL

#** For gcb b cells

samplesToAn = unique(msFNAmerged.goodSamples.gcb@meta.data$sampleName)

for(sm in samplesToAn){
  
  print(sm)
  
  msFNAmerged.goodSamples.gcb.sm = subset(msFNAmerged.goodSamples.gcb, subset = sampleName == sm)
  # first check which clusters in our B cell subclustered dataset are enriched with host lytic markers
  msFNAmerged.goodSamples.gcb.sm.markers <- FindAllMarkers(msFNAmerged.goodSamples.gcb.sm,only.pos = T,logfc.threshold = 0.1)
  msFNAmerged.goodSamples.gcb.sm.markers.sig <- msFNAmerged.goodSamples.gcb.sm.markers[msFNAmerged.goodSamples.gcb.sm.markers$p_val_adj < 0.05,]
  

  totalNumberOfGenes_inLCLsample = nrow(lcl777_b958.ccr)
  numberOflyticGenesAsMarkersInLCLSample = length(host_lyticClusterMarkers)
  
  for(cls in levels(msFNAmerged.goodSamples.gcb.sm.markers.sig$cluster)){
  
  ncells = as.numeric(table(msFNAmerged.goodSamples.gcb.sm$GC_type)[cls]) 
  nCellsInCluster = sum(msFNAmerged.goodSamples.gcb.sm.markers.sig$cluster==cls)
    
  if(nCellsInCluster==0){
    nHostLyticCells_inCluster = 0
    cls_marker_genes = 1
  }else{
    cls_marker_genes = msFNAmerged.goodSamples.gcb.sm.markers.sig$gene[msFNAmerged.goodSamples.gcb.sm.markers.sig$cluster==cls]
    nHostLyticCells_inCluster = length(intersect(host_lyticClusterMarkers,cls_marker_genes))
    
  }
  
  percLytic = nHostLyticCells_inCluster/length(cls_marker_genes)
  
  print(cls)
  print(percLytic)
  print(length(cls_marker_genes))
  
  
  # Fishers exact test to check the proportion of lytic genes in the total set of genes we have in the reference and in our cluster
  
  dat <- data.frame(
    "nonLytic_markers" = c(totalNumberOfGenes_inLCLsample-numberOflyticGenesAsMarkersInLCLSample, length(cls_marker_genes)-nHostLyticCells_inCluster),
    "Lytic_markers" = c(numberOflyticGenesAsMarkersInLCLSample, nHostLyticCells_inCluster),
    row.names = c("LCLsample", "Our_bGC_data_cls"),
    stringsAsFactors = FALSE
  )
  
  ebvLyticmarkers_Enrichment_test <- fisher.test(dat)
  
  print(ebvLyticmarkers_Enrichment_test)
  
  smGrp = as.character(msFNAmerged.goodSamples.gcb.sm$sampleGroup[1])
  totalNumberOfCellsInSample = ncol(msFNAmerged.goodSamples.gcb.sm)
  
  cls_res = c(cls,length(cls_marker_genes),nHostLyticCells_inCluster,percLytic,ebvLyticmarkers_Enrichment_test$p.value,ebvLyticmarkers_Enrichment_test$estimate,paste0(ebvLyticmarkers_Enrichment_test$conf.int,collapse=" - "),sm,totalNumberOfCellsInSample,ncells,smGrp)
  
  msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment.persample = rbind(msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment.persample,cls_res)
  
  
}


}

colnames(msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment.persample) <- c("gcb.cluster","nMarkers","nLyticMarkers","proportionOfLyticMarkers","FishersExact_pvalue","oddsRatio","OddsRatio_Conf.interval","sampleName","totalNumberOfCellsInSample","nCellsInCluster","sampleGroup")


write.csv(msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment.persample,
          file = "msFNAmerged.goodSamples.gcb.clusters.lyticMarkerEnrichment.persample.csv")



#** For Other b cells, here we first prepare a seurat object that contains the plasmablasts along with the B cell types

msFNAmerged.goodSamples = readRDS("data/msFNAmerged.goodSamples.rds")


bcells_insubclustering = colnames(msFNAmerged.goodSamples.memNaivB)
bcells_plasmablasts = rownames(msFNAmerged.goodSamples@meta.data[msFNAmerged.goodSamples@meta.data$rnaclustersvsAdtManAnn2MaxAnn2 == "Plasmablasts",])

intersectingCellRemove = intersect(bcells_insubclustering,bcells_plasmablasts)

bcells_insubclustering = bcells_insubclustering[bcells_insubclustering!=intersectingCellRemove]
bcells_plasmablasts = bcells_plasmablasts[bcells_plasmablasts!=intersectingCellRemove]

selectedBcells = c(bcells_insubclustering,bcells_plasmablasts)

msFNAmerged.goodSamples.memNaivB.plasmablasts = msFNAmerged.goodSamples[,selectedBcells]

sbAnn = sapply(rownames(msFNAmerged.goodSamples.memNaivB.plasmablasts@meta.data),function(rnm){
  an = NA
  if(rnm %in% bcells_insubclustering){
    an = as.character(msFNAmerged.goodSamples.memNaivB@meta.data[rownames(msFNAmerged.goodSamples.memNaivB@meta.data)==rnm,,drop=F]$memnaive_type)
  }else if(rnm %in% bcells_plasmablasts){
    an = as.character(msFNAmerged.goodSamples@meta.data[rownames(msFNAmerged.goodSamples@meta.data)==rnm,,drop=F]$rnaclustersvsAdtManAnn2MaxAnn2)
  }
  an
})

msFNAmerged.goodSamples.memNaivB.plasmablasts@meta.data$subAnnotation = sbAnn

Idents(msFNAmerged.goodSamples.memNaivB.plasmablasts) <- factor(msFNAmerged.goodSamples.memNaivB.plasmablasts@meta.data$subAnnotation)


msFNAmerged.goodSamples.memNaivB.plasmablasts.clusters.lyticMarkerEnrichment.persample = NULL

# Other b cells

samplesToAn2 = unique(msFNAmerged.goodSamples.memNaivB@meta.data$sampleName)

for(sm in samplesToAn2){
  
  print(sm)
  
  msFNAmerged.goodSamples.memNaivB.sm = subset(msFNAmerged.goodSamples.memNaivB.plasmablasts, subset = sampleName == sm)
  # first check which clusters in our B cell subclustered dataset are enriched with host lytic markers
  msFNAmerged.goodSamples.memNaivB.sm.markers <- FindAllMarkers(msFNAmerged.goodSamples.memNaivB.sm,only.pos = T,logfc.threshold = 0.1)
  msFNAmerged.goodSamples.memNaivB.sm.markers.sig <- msFNAmerged.goodSamples.memNaivB.sm.markers[msFNAmerged.goodSamples.memNaivB.sm.markers$p_val_adj < 0.05,]
  
  
  totalNumberOfGenes_inLCLsample = nrow(lcl777_b958.ccr)
  numberOflyticGenesAsMarkersInLCLSample = length(host_lyticClusterMarkers)
  
  for(cls in levels(msFNAmerged.goodSamples.memNaivB.sm.markers.sig$cluster)){
    
    ncells = as.numeric(table(msFNAmerged.goodSamples.memNaivB.sm$subAnnotation)[cls]) 
    nCellsInCluster = sum(msFNAmerged.goodSamples.memNaivB.sm.markers.sig$cluster==cls)
    
    if(nCellsInCluster==0){
      nHostLyticCells_inCluster = 0
      cls_marker_genes = 1
    }else{
      cls_marker_genes = msFNAmerged.goodSamples.memNaivB.sm.markers.sig$gene[msFNAmerged.goodSamples.memNaivB.sm.markers.sig$cluster==cls]
      nHostLyticCells_inCluster = length(intersect(host_lyticClusterMarkers,cls_marker_genes))
      
    }
    
    percLytic = nHostLyticCells_inCluster/length(cls_marker_genes)
    
    print(cls)
    print(percLytic)
    print(length(cls_marker_genes))
    
    
    # Fishers exact test to check the proportion of lytic genes in the total set of genes we have in the reference and in our cluster
    
    dat <- data.frame(
      "nonLytic_markers" = c(totalNumberOfGenes_inLCLsample-numberOflyticGenesAsMarkersInLCLSample, length(cls_marker_genes)-nHostLyticCells_inCluster),
      "Lytic_markers" = c(numberOflyticGenesAsMarkersInLCLSample, nHostLyticCells_inCluster),
      row.names = c("LCLsample", "Our_bGC_data_cls"),
      stringsAsFactors = FALSE
    )
    
    ebvLyticmarkers_Enrichment_test <- fisher.test(dat)
    
    print(ebvLyticmarkers_Enrichment_test)
    
    smGrp = as.character(msFNAmerged.goodSamples.memNaivB.sm$sampleGroup[1])
    totalNumberOfCellsInSample = ncol(msFNAmerged.goodSamples.memNaivB.sm)
    
    cls_res = c(cls,length(cls_marker_genes),nHostLyticCells_inCluster,percLytic,ebvLyticmarkers_Enrichment_test$p.value,ebvLyticmarkers_Enrichment_test$estimate,paste0(ebvLyticmarkers_Enrichment_test$conf.int,collapse=" - "),sm,totalNumberOfCellsInSample,ncells,smGrp)
    
    msFNAmerged.goodSamples.memNaivB.plasmablasts.clusters.lyticMarkerEnrichment.persample = rbind(msFNAmerged.goodSamples.memNaivB.plasmablasts.clusters.lyticMarkerEnrichment.persample,cls_res)
    
    
  }
  
  
}

colnames(msFNAmerged.goodSamples.memNaivB.plasmablasts.clusters.lyticMarkerEnrichment.persample) <- c("bcell.cluster","nMarkers","nLyticMarkers","proportionOfLyticMarkers","FishersExact_pvalue","oddsRatio","OddsRatio_Conf.interval","sampleName","totalNumberOfCellsInSample","nCellsInCluster","sampleGroup")



write.csv(msFNAmerged.goodSamples.memNaivB.plasmablasts.clusters.lyticMarkerEnrichment.persample,
          file = "msFNAmerged.goodSamples.memNaivB.plasmablasts.clusters.lyticMarkerEnrichment.persample.csv")



# plot host lytic gene enrichment for memNaivB.plasmablasts (for gc bcells numbers are too low)

memNaivB.pb.clusters.lyticMarkerEnrichment = as.data.frame(msFNAmerged.goodSamples.memNaivB.plasmablasts.clusters.lyticMarkerEnrichment.persample)
memNaivB.pb.clusters.lyticMarkerEnrichment$Bcellstype = "Extrafollicular"
colnames(memNaivB.pb.clusters.lyticMarkerEnrichment)[1] <- "cluster"

memNaivB.pb.clusters.lyticMarkerEnrichment$oddsRatio = as.numeric(memNaivB.pb.clusters.lyticMarkerEnrichment$oddsRatio)
memNaivB.pb.clusters.lyticMarkerEnrichment$p.value = as.numeric(memNaivB.pb.clusters.lyticMarkerEnrichment$FishersExact_pvalue) < 0.05

ggplot(data = memNaivB.pb.clusters.lyticMarkerEnrichment, aes(x = cluster, y = oddsRatio,fill=p.value)) + coord_flip() + 
  geom_point(shape = 21, colour = "black",size=4) + 
  geom_hline(yintercept=1, linetype="dashed", color = "grey") + 
  scale_fill_manual(values = c("white","red")) +
  theme_classic() + 
  ylab("OR") + guides(fill=guide_legend(title="p.value < 0.05")) +
  facet_grid(cols  = vars(sampleName))

ggsave("Extrafollicular_memNaivB.pb.clusters.lyticMarker_overrepresentation.pdf",
       height = 100, width = 400, units = "mm")




## mapping of the public LCL cells to our B clusters, and identify which cluster is enriched with cells with LCL like transcriptomes ##

# B clusters

bcellClusters.lcl777_b958.anchors <- FindTransferAnchors(
  reference = msFNAmerged.goodSamples.memNaivB,
  query = lcl777_b958.ccr,
  dims = 1:50
)

bcellClusters.predictions <- TransferData(anchorset = bcellClusters.lcl777_b958.anchors, 
                                          refdata = msFNAmerged.goodSamples.memNaivB$memnaive_type)



lcl777_b958.ccr.bcluster.predFreq = table(bcellClusters.predictions$predicted.id)/sum(table(bcellClusters.predictions$predicted.id))


lcl777_b958.ccr[["bcellClusters.id"]] <-  bcellClusters.predictions$predicted.id
lcl777_b958.ccr[["bcellClusters.scoreMax"]] <-  bcellClusters.predictions$prediction.score.max


lcl777_clusters_prediction_prop = table(lcl777_b958.ccr$seurat_clusters,lcl777_b958.ccr$bcellClusters.id)
lcl777_clusters_prediction_prop = t(apply(lcl777_clusters_prediction_prop,1,function(x) 100 * (x/sum(x))))

pdf("lcl777_bcellclusters_prediction_prop_heatmap.pdf")
pheatmap::pheatmap(lcl777_clusters_prediction_prop,scale="none",cluster_cols = F)
dev.off()

bcluster.Freq = table(msFNAmerged.goodSamples.memNaivB$memnaive_type) / sum(table(msFNAmerged.goodSamples.memNaivB$memnaive_type))

# We mapped LCL cells from lcl777 their closest b cell subset in our dataset using our data as reference.  26% of the 1643
# lcl cells in lcl 777 were predicted to be closest to the DN population (which made up just 2.6% of our b cells population). Interestingly,
# the more lytic lcl clusters in lcl 777 were overwhelmingly predicted to be similar transcription ally to DN, for instance 83% cells (30 of 36) of cluster 7 of lcl 777 were predicted to be closest to DN. Additionally,
# DN predictions from lcl 777 cluster 7 had highest median mapping prediction scores compared to DN prediction from other clusters of lcl 777, 
# suggesting overall that the DN cluster is the closest transcriptionally to LCLs, and particularly lytic lcls.

ggplot(lcl777_b958.ccr@meta.data, aes(x=factor(bcellClusters.id), y=bcellClusters.scoreMax,color = factor(seurat_clusters)))+
  geom_boxplot() +  theme(legend.position = "top") +  theme_minimal() + xlab("") + ylab("Mapping prediction Score")

ggsave("lcl777_bcellclusters_prediction_scoreDistributions_bylcl777cluster.pdf",height = 100, width = 200, units = "mm")

bcellClusters.predictionScoresByCellTypeAndCluster <- lcl777_b958.ccr@meta.data %>% group_by(bcellClusters.id,seurat_clusters) %>% 
  summarise(n=n(),meanPredictionScore = mean(bcellClusters.scoreMax))


# map the second lcl data lcl777_m81.ccr to B clusters 

bcellClusters.lcl777_m81.ccr.anchors <- FindTransferAnchors(
  reference = msFNAmerged.goodSamples.memNaivB,
  query = lcl777_m81.ccr,
  dims = 1:50
)

lcl777_m81.ccr.bcellClusters.predictions <- TransferData(anchorset = bcellClusters.lcl777_m81.ccr.anchors, 
                                          refdata = msFNAmerged.goodSamples.memNaivB$memnaive_type)



lcl777_m81.ccr.bcluster.predFreq = table(lcl777_m81.ccr.bcellClusters.predictions$predicted.id)/sum(table(lcl777_m81.ccr.bcellClusters.predictions$predicted.id))


## map the third lcl sample (lcl461_b958.ccr)

bcellClusters.lcl461_b958.ccr.anchors <- FindTransferAnchors(
  reference = msFNAmerged.goodSamples.memNaivB,
  query = lcl461_b958.renamed.ccr,
  dims = 1:50
)

lcl461_b958.bcellClusters.predictions <- TransferData(anchorset = bcellClusters.lcl461_b958.ccr.anchors, 
                                                         refdata = msFNAmerged.goodSamples.memNaivB$memnaive_type,k.weight=30)



lcl461_b958.ccr.bcluster.predFreq = table(lcl461_b958.bcellClusters.predictions$predicted.id)/sum(table(lcl461_b958.bcellClusters.predictions$predicted.id))


## map the fourth lcl sample (gm12878.ccr, public data)

bcellClusters.gm12878.ccr.anchors <- FindTransferAnchors(
  reference = msFNAmerged.goodSamples.memNaivB,
  query = gm12878.ccr,
  dims = 1:50
)

gm12878.ccr.bcellClusters.predictions <- TransferData(anchorset = bcellClusters.gm12878.ccr.anchors, 
                                                      refdata = msFNAmerged.goodSamples.memNaivB$memnaive_type)



gm12878.ccr.bcluster.predFreq = table(gm12878.ccr.bcellClusters.predictions$predicted.id)/sum(table(gm12878.ccr.bcellClusters.predictions$predicted.id))


## map the fifth lcl sample (gm12878.ccr, public data)

bcellClusters.gm18502.ccr.anchors <- FindTransferAnchors(
  reference = msFNAmerged.goodSamples.memNaivB,
  query = gm18502.ccr,
  dims = 1:50
)

gm18502.ccr.bcellClusters.predictions <- TransferData(anchorset = bcellClusters.gm18502.ccr.anchors, 
                                                      refdata = msFNAmerged.goodSamples.memNaivB$memnaive_type)



gm18502.ccr.bcluster.predFreq = table(gm18502.ccr.bcellClusters.predictions$predicted.id)/sum(table(gm18502.ccr.bcellClusters.predictions$predicted.id))


## Conclusion: DN population of B cells looks the most similar to LCL cells from all lcl samples.
## and this is the case in all 5 LCL samples from : https://elifesciences.org/articles/62586#s2

bcluster.Freq.allSamples = data.frame(bcluster.Freq)
colnames(bcluster.Freq.allSamples) <- c("Bcell_type","pCurrentData")

bcluster.Freq.allSamples$pLCL777_b958 <- 0
bcluster.Freq.allSamples$pLCL777_b958[match(names(lcl777_b958.ccr.bcluster.predFreq),bcluster.Freq.allSamples$Bcell_type)] <- lcl777_b958.ccr.bcluster.predFreq

bcluster.Freq.allSamples$pLCL777_m81 <- 0
bcluster.Freq.allSamples$pLCL777_m81[match(names(lcl777_m81.ccr.bcluster.predFreq),bcluster.Freq.allSamples$Bcell_type)] <- lcl777_m81.ccr.bcluster.predFreq

bcluster.Freq.allSamples$pLCL461_b958 <- 0
bcluster.Freq.allSamples$pLCL461_b958[match(names(lcl461_b958.ccr.bcluster.predFreq),bcluster.Freq.allSamples$Bcell_type)] <- lcl461_b958.ccr.bcluster.predFreq

bcluster.Freq.allSamples$pLCL_gm12878 <- 0
bcluster.Freq.allSamples$pLCL_gm12878[match(names(gm12878.ccr.bcluster.predFreq),bcluster.Freq.allSamples$Bcell_type)] <- gm12878.ccr.bcluster.predFreq

bcluster.Freq.allSamples$pLCL_gm18502 <- 0
bcluster.Freq.allSamples$pLCL_gm18502[match(names(gm18502.ccr.bcluster.predFreq),bcluster.Freq.allSamples$Bcell_type)] <- gm18502.ccr.bcluster.predFreq

rownames(bcluster.Freq.allSamples) <- bcluster.Freq.allSamples$Bcell_type

bcluster.Freq.allSamples.prop <- bcluster.Freq.allSamples[,-1] * 100

pdf("bcluster_prediction_Freq_allSamples_prop_heatmap.pdf")
pheatmap::pheatmap(bcluster.Freq.allSamples.prop,scale="none",cluster_cols = F,cluster_rows = T)
dev.off()



# gc

bgcCls.lcl777_b958.anchors <- FindTransferAnchors(
  reference = msFNAmerged.goodSamples.gcb,
  query = lcl777_b958.ccr,
  dims = 1:50
)

bgcCls.predictions <- TransferData(anchorset = bgcCls.lcl777_b958.anchors, 
                                          refdata = msFNAmerged.goodSamples.gcb$GC_type)



lcl777_b958.ccr[["bgcClusters.id"]] <-  bgcCls.predictions$predicted.id
lcl777_b958.ccr[["bgcClusters.scoreMax"]] <-  bgcCls.predictions$prediction.score.max


barplot(table(bgcCls.predictions$predicted.id)/sum(table(bgcCls.predictions$predicted.id)),ylim=c(0,0.6))


lcl777_bgcclusters_prediction_prop = table(lcl777_b958.ccr$seurat_clusters,lcl777_b958.ccr$bgcClusters.id)
lcl777_bgcclusters_prediction_prop = t(apply(lcl777_bgcclusters_prediction_prop,1,function(x) 100 * (x/sum(x))))

pdf("lcl777_bgcclusters_prediction_prop_heatmap.pdf")
pheatmap::pheatmap(lcl777_bgcclusters_prediction_prop,scale="none",cluster_cols = F)
dev.off()

# Mapping of the same cells to our gc bcell data showed majority were predicted to be closest to the prolif.T cluster. 
# According to this paper (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7469576/), ebv transformed b cells significantly upregulate CD3D (compared to normal B cells)
# Another paper also shows B cells with T cell like features in ebv infection situation (in a study of EBV-Associated Gastric Cancers)
# This b cell cluster with t cell features expressed CD3D compared to other b cell clusters (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10296402/figure/cancers-15-03178-f004/)
# this prolif.T just makes up 15 percent of our gc b cells, but 71% of the lcl 777 cells were close to this cluster transcriptionally.
# thus this gc b cell cluster with T cell like features is likely to be an ebv-transformed b cell population.

# the prediction of prolif.T overall had a higher prediction score compared to lcl cells mapped to be GC.B and LZ types, suggesting that
# within the gc b cell subset, the transcriptional similarity of LCL cells is overall closest to the prolif.T cells in our gc b cells, and then to LZ and GC.B, but not to any dark zone cells in our data.

table(msFNAmerged.goodSamples.gcb$GC_type) / sum(table(msFNAmerged.goodSamples.gcb$GC_type))


ggplot(lcl777_b958.ccr@meta.data, aes(x=factor(bgcClusters.id), y=bgcClusters.scoreMax,color = factor(seurat_clusters)))+
  geom_boxplot() +  theme(legend.position = "top") +  theme_minimal() + xlab("") + ylab("Mapping prediction Score")

ggsave("lcl777_gcbcellclusters_prediction_scoreDistributions_bylcl777cluster.pdf",height = 100, width = 200, units = "mm")



# which of the bcell clusters do gc b cells resemble most transcriptionally ?? 
# map gc cells to the bclusters

gc_to_bclust.anchors <- FindTransferAnchors(
  reference = msFNAmerged.goodSamples.memNaivB,
  query = msFNAmerged.goodSamples.gcb,
  dims = 1:50
)

gc_to_bclust.predictions <- TransferData(anchorset = gc_to_bclust.anchors, 
                                   refdata = msFNAmerged.goodSamples.memNaivB$memnaive_type)


gc_to_bclust.predFreq = table(gc_to_bclust.predictions$predicted.id)/sum(table(gc_to_bclust.predictions$predicted.id))


msFNAmerged.goodSamples.gcb.bclustPred = msFNAmerged.goodSamples.gcb

msFNAmerged.goodSamples.gcb.bclustPred[["bcellClusters.id"]] <-  gc_to_bclust.predictions$predicted.id
msFNAmerged.goodSamples.gcb.bclustPred[["bcellClusters.scoreMax"]] <-  gc_to_bclust.predictions$prediction.score.max


msFNAmerged.goodSamples.gcb.bclust_prediction = table(msFNAmerged.goodSamples.gcb.bclustPred$GC_type,msFNAmerged.goodSamples.gcb.bclustPred$bcellClusters.id)
msFNAmerged.goodSamples.gcb.bclust_prediction_prop = t(apply(msFNAmerged.goodSamples.gcb.bclust_prediction,1,function(x) 100 * (x/sum(x))))

pdf("gc_to_bclust_prediction_prop_heatmap.pdf")
pheatmap::pheatmap(msFNAmerged.goodSamples.gcb.bclust_prediction_prop,scale="none",cluster_cols = F)
dev.off()


# we map gc bcells to other b cells including plasmablasts
# and no profil.T
#msFNAmerged.goodSamples.memNaivB.plasmablasts@meta.data$subAnnotation

msFNAmerged.goodSamples.gcb.noProfilT = subset(msFNAmerged.goodSamples.gcb, subset= GC_type %in% c("LZ","DZ.1","DZ.2","GC.B"))

gc_to_bclust.anchors.2 <- FindTransferAnchors(
  reference = msFNAmerged.goodSamples.memNaivB.plasmablasts,
  query = msFNAmerged.goodSamples.gcb.noProfilT,
  dims = 1:50
)

gc_to_bclust.predictions.2 <- TransferData(anchorset = gc_to_bclust.anchors.2, 
                                         refdata = msFNAmerged.goodSamples.memNaivB.plasmablasts$subAnnotation)


gc_to_bclust.predFreq.2 = table(gc_to_bclust.predictions.2$predicted.id)/sum(table(gc_to_bclust.predictions.2$predicted.id))


msFNAmerged.goodSamples.gcb.bclustPred2 = msFNAmerged.goodSamples.gcb.noProfilT

msFNAmerged.goodSamples.gcb.bclustPred2[["bcellClusters.id"]] <-  gc_to_bclust.predictions.2$predicted.id
msFNAmerged.goodSamples.gcb.bclustPred2[["bcellClusters.scoreMax"]] <-  gc_to_bclust.predictions.2$prediction.score.max


msFNAmerged.goodSamples.gcb.bclust_prediction2 = table(msFNAmerged.goodSamples.gcb.bclustPred2$GC_type,msFNAmerged.goodSamples.gcb.bclustPred2$bcellClusters.id)
msFNAmerged.goodSamples.gcb.bclust_prediction_prop2 = t(apply(msFNAmerged.goodSamples.gcb.bclust_prediction2[-5,],1,function(x) 100 * (x/sum(x))))

pdf("gc_to_bclust_plasmablasts_prediction_prop_heatmap.pdf")
pheatmap::pheatmap(msFNAmerged.goodSamples.gcb.bclust_prediction_prop2,scale="none",cluster_cols = F)
dev.off()


write.csv(msFNAmerged.goodSamples.gcb.bclust_prediction_prop2,
          file = "gc_to_bclust_plasmablasts_prediction_prop_heatmap.csv")


# calculate the frequencies by status

msFNAmerged.goodSamples.gcb.bclustPred2.ms = subset(msFNAmerged.goodSamples.gcb.bclustPred2,subset=sampleGroup == "MS")
msFNAmerged.goodSamples.gcb.bclust_prediction2.ms = table(msFNAmerged.goodSamples.gcb.bclustPred2.ms$GC_type,msFNAmerged.goodSamples.gcb.bclustPred2.ms$bcellClusters.id)
msFNAmerged.goodSamples.gcb.bclust_prediction_prop2.ms = t(apply(msFNAmerged.goodSamples.gcb.bclust_prediction2.ms[-5,],1,function(x) 100 * (x/sum(x))))


msFNAmerged.goodSamples.gcb.bclustPred2.ctrl = subset(msFNAmerged.goodSamples.gcb.bclustPred2,subset=sampleGroup == "CTRL")
msFNAmerged.goodSamples.gcb.bclust_prediction2.ctrl = table(msFNAmerged.goodSamples.gcb.bclustPred2.ctrl$GC_type,msFNAmerged.goodSamples.gcb.bclustPred2.ctrl$bcellClusters.id)
msFNAmerged.goodSamples.gcb.bclust_prediction_prop2.ctrl = t(apply(msFNAmerged.goodSamples.gcb.bclust_prediction2.ctrl[-5,],1,function(x) 100 * (x/sum(x))))







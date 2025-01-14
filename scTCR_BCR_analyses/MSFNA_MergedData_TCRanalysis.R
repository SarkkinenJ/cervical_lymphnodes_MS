library(Seurat)
library(ggplot2)
library(plyr)
library(ggpubr)
library(dplyr)


library(scRepertoire)
library("reticulate")

source("TCR_functions.R")


# load the merged Seurat object of the entire data ####

msFNAmerged <- readRDS(file = "data/msFNAmerged.rds")


# read the filtered TCR files into a list for TCR analysis #### 


# Use the directory where the sc TCR files are found for all samples (example naming: MS002_filtered_contig_annotations.csv)

filteredContigFileNames = list.files("/scRNAseq/vdjdata/",pattern ="_contig_",full.names=T,recursive = TRUE)



MS5_TCRcontigList <- lapply(filteredContigFileNames, function(x) read.csv(file = x,header=T))
names(MS5_TCRcontigList) <- sapply(filteredContigFileNames,
                                   function(y) paste0(unlist(strsplit(basename(y),"_"))[1],collapse="_"))



# add sample barcode prefixes to the tcr data 
for(x in 1:length(MS5_TCRcontigList)){
 MS5_TCRcontigList[[x]]$barcode <- paste(names(MS5_TCRcontigList[x]),MS5_TCRcontigList[[x]]$barcode,sep="_")
}



# Combine the Contigs (of alpha and beta chain sequences for each cell are combined)

sampleGrp <- msFNAmerged@meta.data$sampleGroup[match(names(MS5_TCRcontigList),msFNAmerged@meta.data$orig.ident)]
sampleNms <- msFNAmerged@meta.data$sampleName[match(names(MS5_TCRcontigList),msFNAmerged@meta.data$orig.ident)]


combined_MS5_TCRcontigList <- combineTCR(MS5_TCRcontigList, 
                                         samples = sampleNms, 
                                         ID = sampleGrp, cells ="T-AB")


# modify the barcodes to match the barcodes in the seurat object 
for(x in 1:length(combined_MS5_TCRcontigList)){
  combined_MS5_TCRcontigList[[x]]$barcode <- sapply(combined_MS5_TCRcontigList[[x]]$barcode,
                                                    function(x) paste0(unlist(strsplit(x,"_"))[c(-2,-3)],collapse="_"))
}
names(combined_MS5_TCRcontigList) <- sampleNms



# combine to seurat data 
# Samples "MS001","CTRL012","CTRL013","CTRL014" are removed from analyses because they contained significantly less number of cells compared to other samples

samnamesSelected = samnames[!samnames %in% c("MS001","CTRL012","CTRL013","CTRL014")]
msFNAmerged.fortcr <- subset(x = msFNAmerged, subset = sampleName %in% samnamesSelected)

msFNA.merged.tcr <- combineExpression(combined_MS5_TCRcontigList, msFNAmerged.fortcr, 
                                      cloneCall="gene+nt", groupBy = "sample")




#### Save important data objects ####
saveRDS(combined_MS5_TCRcontigList, file = "data/combined_MS5_TCRcontigList.rds")
saveRDS(msFNA.merged.tcr, file = "data/msFNA.merged.tcr.rds")




### Complete TCR diversity analysis for all T cell cell types ####

tcrCells = rownames(msFNA.merged.tcr@meta.data[!is.na(msFNA.merged.tcr@meta.data$barcode),])

msFNA.merged.tcrData = subset(msFNA.merged.tcr, subset = barcode %in% tcrCells)

cellTypesWithTCRs <- names(table(msFNA.merged.tcrData@meta.data$rnaclustersvsAdtManAnn2MaxAnn2))
ResultDir = "TCRanalysis_Results_plotsWithSampleNames/"

for(ctcr in c("allTCR",cellTypesWithTCRs)){
  
  print(ctcr)
  if(ctcr == "allTCR"){
    msFNA.merged.bcrTempData = msFNA.merged.tcrData
  }else{
    msFNA.merged.tcrTempData = subset(msFNA.merged.tcrData,subset = rnaclustersvsAdtManAnn2MaxAnn2 == ctcr)
  }
  
  #msFNA.merged.tcrTempData = subset(msFNA.merged.tcrData,subset = rnaclustersvsAdtManAnn2MaxAnn2 == ctcr)
  combined_MS5_tcrcontigList_tcrTempData <- expression2List(msFNA.merged.tcrTempData,group="sampleName")
  
  write.to = paste0(ResultDir,ctcr,"_finalAnnotation/")
  dir.create(write.to)
  
  
  # # % of unique clonotypes 
  # quantContig(combined_MS5_tcrcontigList_tcrTempData, cloneCall="CTstrict", scale = T)
  # ggsave(paste0(write.to,ctcr,"_percentOfUniqueClonotypes.pdf"), height = 210, width = 320, units = "mm")
  # 
  # # % of unique clonotypes by group
  # quantContig(combined_MS5_tcrcontigList_tcrTempData, cloneCall="CTstrict", scale = T,group="sampleGroup")
  # ggsave(paste0(write.to,ctcr,"_percentOfUniqueClonotypes_ByGroup.pdf"), height = 210, width = 320, units = "mm")
  # 
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
    g + geom_boxplot(outlier.colour = NA) +  geom_text(label=percOfuniqd$values,size = 1,colour = "red",check_overlap = TRUE, position=position_jitter(width=0.50)) +
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
    g + geom_boxplot(outlier.colour = NA) +  geom_text(label=percOfuniqd$values,size = 1,colour = "red",check_overlap = TRUE, position=position_jitter(width=0.50)) + 
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
      g + geom_boxplot(outlier.colour = NA) +  geom_text(label=diversitytcrClones$sample,size = 1,colour = "red",check_overlap = TRUE, position=position_jitter(width=0.50)) +
        #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
        #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
        #scale_color_manual(values=c("lightblue","Brown")) + 
        scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
        geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
        theme_classic() + theme(legend.position = "none") + 
        theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
        ylab("Shannon diversity") + ggtitle(ctcr) +  NoLegend() + 
      
      #stat_pvalue_manual(stat.test, x = "timePoint", hide.ns = TRUE,label = "p.signif",label.size=6,inherit.aes = FALSE,tip.length = 0)
      
      # significance with pvalues
      stat_pvalue_manual(stat.test, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)
      ggsave(paste0(write.to,ctcr,"_tcrdiversity_Shannon.pdf"), width = 2, height = 4)
      
  
    }else{

      g <- ggplot(data = diversitytcrClones, aes(x = sampleGroup,y = Shannon,fill=sampleGroup))  
      g + geom_boxplot(outlier.colour = NA) +  geom_text(label=diversitytcrClones$sample,size = 1,colour = "red",check_overlap = TRUE, position=position_jitter(width=0.50)) +
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
  diversityValuesAAtemp$sample <- c(names(repsTemp[grepl("MS",names(repsTemp))]),names(repsTemp[grepl("CTRL",names(repsTemp))]))
  
  diversityValuesAAtemp_statp  = compare_means(diversity ~ sampleGroup, diversityValuesAAtemp,method="t.test",p.adjust.method="BH")
  
  stat.test <- diversityValuesAAtemp_statp %>%
    mutate(y.position = rep(max(diversityValuesAAtemp$diversity) + 0.01,1))
  
  
  g <- ggplot(data = diversityValuesAAtemp, aes(x = sampleGroup,y = diversity,fill=sampleGroup))  
  g + geom_boxplot(outlier.colour = NA) +  geom_text(label=diversityValuesAAtemp$sample,size = 1,colour = "red",check_overlap = TRUE, position=position_jitter(width=0.50)) +
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




## work with database scan results ####

## First scan the VDJdb  TCR database for similar sequences using runTCRdbScan.sh (which runs tcr_VDJDBScan.R for all samples at once)
## the results from scanning VDJdb are directory named: dbScanResults

msFNA.merged.tcr@meta.data$predictedSpecificityVDJdbP <- NA
msFNA.merged.tcr@meta.data$predictedEpitopeGeneVDJdbP <- NA

# add specificity after hla matching of the cell type with the hla of the tcr on the database
msFNA.merged.tcr@meta.data$predictedSpecificity_hlamatched_VDJdbP <- NA
msFNA.merged.tcr@meta.data$predictedEpitopeGene_hlamatched_VDJdbP <- NA
msFNA.merged.tcr@meta.data$queryCelltype_targetTCRHLA <- NA


numberOfUniqueClonesWithEBVhits <- c()
numberOfUniqueClonesWithEBVhitsPer1000 <- c()
names(numberOfTCRs) <- names(combined_MS5_TCRcontigList)
EnrichmentOfCloneType <- NULL
nHitsPerPathType <- NULL
EBVTargettingCellTypes <- NULL
nCells <- length(unique(msFNA.merged.tcr@meta.data$rnaclustersvsAdtManAnn2MaxAnn2))
cnames <- unique(msFNA.merged.tcr@meta.data$rnaclustersvsAdtManAnn2MaxAnn2)

scannedDB="_paired_VDJdb_scan.csv"
#scannedDB="_beta_VDJdb_scan.csv"
#ebvPath="Epstein Barr virus (EBV)"
ebvPath="EBV"
#pathCol= "Pathology"
pathCol="antigen.species"

for(sam in names(combined_MS5_TCRcontigList)){
  samRes = read.csv2(paste0("dbScanResults/",sam,scannedDB),sep=",")
  samRes <- samRes[!is.na(samRes$antigen.species),]
  
  bestHitData = NULL
  exactHitData = NULL
  #samRes <- samRes[samRes$combinedPaired_tcrdist <= 24,]
  bestHitDatatemp <- sapply(unique(samRes$barcode),function(x){
    
    hitsForBarcode = samRes[samRes$barcode==x,]
    
    # test if HLA is for the hit TCR is the expected HLA (class 1 for cd8 cells and class 2 for other CD4s); added after Jani's comments
    # that tcr hits should have the expected HLA presentation, that is class 1 presented epitope associated TCRs should be CD8s
    
    cellt <- msFNA.merged.tcr@meta.data[match(x,msFNA.merged.tcr@meta.data$barcode),]$rnaclustersvsAdtManAnn2MaxAnn2
    cellt_hla = ifelse(cellt %in% c("Memory CD8","Naive CD8"),"MHCI","MHCII")
    
    hitsForBarcode$queryCelltype_targetTCRHLA = paste(cellt,hitsForBarcode$mhc.class,sep="_")
    
    hitsForBarcode = hitsForBarcode[hitsForBarcode$mhc.class==cellt_hla,]
    exacth = hitsForBarcode[hitsForBarcode$combinedPaired_tcrdist==0,][1,]
    
    
    hitsForBarcode = hitsForBarcode[order(hitsForBarcode$combinedPaired_tcrdist,hitsForBarcode$cdr3b_tcrdist,hitsForBarcode$cdr3a_tcrdist),] # order them by combined dist and then b dist and a dist, take the top one
    
    bestHitData <<- rbind(bestHitData,hitsForBarcode[1,])
    
    exactHitData <<- rbind(exactHitData,exacth)
    
    hitsForBarcode[1,]
  })
  
  bestHitData <- bestHitData[!is.na(bestHitData$antigen.species),]
  exactHitData <- exactHitData[!is.na(exactHitData$antigen.species),]
  
  
  write.csv(bestHitData,paste0("dbScanResults/",sam,"_bestHitsPerClone_hlaMatched",scannedDB))
  write.csv(exactHitData,paste0("dbScanResults/",sam,"_exactHitsPerClone_hlaMatched",scannedDB))
  
  samRes = bestHitData
  
  scoreDistributionForPath = table(samRes$combinedPaired_tcrdist,samRes$antigen.species)
  write.csv(scoreDistributionForPath,file=paste0("dbScanResults/",sam,"_VDJpairedScan_bestHits__hlaMatched_PathologyScoreDistribution.csv"))
  
  samResEBV <- samRes[samRes[,pathCol]==ebvPath,]
  uniqClonesForEBV <- length(unique(samResEBV$barcode))
  numberOfUniqueClonesWithEBVhits <- c(numberOfUniqueClonesWithEBVhits,uniqClonesForEBV)
  
  numsampTCRs <- nrow(combined_MS5_TCRcontigList[[sam]])
  
  numberOfUniqueClonesWithEBVhitsPer1000 <- c(numberOfUniqueClonesWithEBVhitsPer1000,(uniqClonesForEBV/numsampTCRs) * 1000)
  print(sam)
  
  # add predicted specificity to the seurat meta data of the cells.
  sampleTCRsNotFilteredOutBySeurat <- samRes[samRes$barcode %in% rownames(msFNA.merged.tcr@meta.data),]
  msFNA.merged.tcr@meta.data[match(sampleTCRsNotFilteredOutBySeurat$barcode,rownames(msFNA.merged.tcr@meta.data)),]$predictedSpecificity_hlamatched_VDJdbP <- sampleTCRsNotFilteredOutBySeurat$antigen.species
  msFNA.merged.tcr@meta.data[match(sampleTCRsNotFilteredOutBySeurat$barcode,rownames(msFNA.merged.tcr@meta.data)),]$predictedEpitopeGene_hlamatched_VDJdbP <- sampleTCRsNotFilteredOutBySeurat$antigen.gene
  
  msFNA.merged.tcr@meta.data[match(sampleTCRsNotFilteredOutBySeurat$barcode,rownames(msFNA.merged.tcr@meta.data)),]$queryCelltype_targetTCRHLA <- sampleTCRsNotFilteredOutBySeurat$queryCelltype_targetTCRHLA
  
  
  # ebv specificity rates per cell type
  cellTEBV <- table(msFNA.merged.tcr@meta.data[unique(samResEBV$barcode),]$rnaclustersvsAdtManAnn2MaxAnn2)/uniqClonesForEBV
  EBVTargettingCellTypes <- cbind(EBVTargettingCellTypes,sapply(unique(msFNA.merged.tcr@meta.data$rnaclustersvsAdtManAnn2MaxAnn2),function(x) ifelse(x %in% names(cellTEBV),cellTEBV[x],0)))
  
  nTCRsPerSamp <- c()
  # nums to compare For all hit types
  for(pat in unique(VDJdblatest$antigen.species)){
    
    samResP <- samRes[samRes$antigen.species==pat,]
    
    numberOfEBVhits <- ifelse(nrow(samResP) > 0, length(unique(samResP$barcode)), 0)
    numberOfNonEBVhits <- length(unique(samRes[samRes$antigen.species!=pat,]$barcode))
    
    # ebv in database
    numberOfEBVrecordsinDB <- sum(VDJdblatest$antigen.species==pat,na.rm=T)
    numberOfNonEBVrecordsinDB <- sum(VDJdblatest$antigen.species!=pat,na.rm=T)
    
    dForF <- data.frame(hits = c(numberOfEBVhits,numberOfEBVrecordsinDB), nonhits = c(numberOfNonEBVhits,numberOfNonEBVrecordsinDB))
    rownames(dForF) <- c("InSample","InReferenceDB")
    
    fres = fisher.test(dForF)
    
    
    EnrichmentOfCloneType = rbind(EnrichmentOfCloneType,c(sam,pat,fres$p.value,fres$estimate,paste(fres$conf.int,collapse="-")))
    
    nTCRsPerSamp <- c(nTCRsPerSamp,(numberOfEBVhits/numsampTCRs) * 1000)
    
  }
  
  nHitsPerPathType <- cbind(nHitsPerPathType,nTCRsPerSamp)
  
}

rownames(EBVTargettingCellTypes) <- cnames
colnames(EBVTargettingCellTypes) <- names(combined_MS5_TCRcontigList)
EBVTargettingCellTypesComp = t(apply(EBVTargettingCellTypes,1,function(x) c(wilcox.test(x[c(1,2,3,4,6,7)],x[c(5,8,9)])$p.value,
                                                                        log2(mean(x[c(1,2,3,4,6,7)])/mean(x[c(5,8,9)])))))

colnames(EBVTargettingCellTypesComp) <- c("unpaired two-samples wilcox.test.p.value","log2FC(meanMS/meanCTRL)")

EBVTargettingCellTypes <- cbind(EBVTargettingCellTypes,EBVTargettingCellTypesComp)


write.csv(EBVTargettingCellTypes,file="dbScanResults/VDJpairedScan_besthits_hlaMatched_EBVTargettingCellTypes.csv")


rownames(nHitsPerPathType) <- unique(VDJdblatest$antigen.species)
colnames(nHitsPerPathType) <- names(combined_MS5_TCRcontigList)
nHitsPerPathTypestatComp = t(apply(nHitsPerPathType + 1,1,function(x) c(wilcox.test(x[c(1,2,3,4,6,7)],x[c(5,8,9)])$p.value,
                                                                        log2(mean(x[c(1,2,3,4,6,7)])/mean(x[c(5,8,9)])))))


colnames(nHitsPerPathTypestatComp) <- c("unpaired two-samples wilcox.test.p.value","log2FC(meanMS/meanCTRL)")

nHitsPerPathTypeFull = cbind(nHitsPerPathType,nHitsPerPathTypestatComp)
write.csv(nHitsPerPathTypeFull,file="dbScanResults/VDJpairedScan_besthits_hlaMatched_nHitsPerPathTypeFull.csv")


colnames(EnrichmentOfCloneType) <- c("sample","pathology","pvalue","OddsRatio","OR95ConfInt")
padj = c()
for (sam in names(combined_MS5_TCRcontigList)){
  pvals = EnrichmentOfCloneType[EnrichmentOfCloneType[,1]==sam,][,3]
  pa = round(p.adjust(as.numeric(pvals), "BH"), 3)
  padj = c(padj,pa)
}

EnrichmentOfCloneType <- cbind(EnrichmentOfCloneType,padj)


# pathologies ordered by number of hits
meanHitForPath = apply(nHitsPerPathTypeFull[,1:9],1,mean)

meanHitForPath_ordered = sort(meanHitForPath,decreasing=T)
meanHitForPath_ordered_besthits = sort(meanHitForPath,decreasing=T)

barplot(meanHitForPath_ordered/sum(meanHitForPath_ordered), cex.names=0.8,las=2)

rank(meanHitForPath_ordered)

numberOfPathRecordsInDb_ordered = sort(table(VDJdblatest$antigen.species),decreasing=T)
barplot(numberOfPathRecordsInDb_ordered/sum(numberOfPathRecordsInDb_ordered), cex.names=0.8,las=2)

rank(numberOfPathRecordsInDb_ordered)

# combined barplot
numPathIndb = table(VDJdblatest$antigen.species)/sum(table(VDJdblatest$antigen.species))

numPathsInHit = meanHitForPath/sum(meanHitForPath)
numPathsInHit = numPathsInHit[names(numPathIndb)]

numPathCombined = rbind(numPathIndb,numPathsInHit)

barplot(numPathCombined, col=c("darkblue","red"), ylim=c(0,0.8),
        legend = c("Proportion in VDJdb","Mean Proportion of predicited specificity (all samples)"), beside=TRUE,cex.names=0.8,las=2)

numPathsRankRatio = rank(numPathsInHit)/rank(numPathIndb)
numPathsRankRatio[which(numPathsInHit > 0)][which.max(numPathsRankRatio[which(numPathsInHit > 0)])]


numPathCombinedRank = rbind(rank(numPathIndb),rank(numPathsInHit))

barplot(numPathCombinedRank, col=c("darkblue","red"),ylim=c(1,30),
        legend = c("Proportion in VDJdb","Mean Proportion of predicited specificity (all samples)"), beside=TRUE,cex.names=0.8,las=2)





## Evaluate cells by expansion types ####

# Definition for expansion types (https://jitc.bmj.com/content/10/6/e004512)
# Hyperexpanded (when x>0.01), large (0.001 < x <= 0.01), medium (10^-4 < x <= 0.001), and small/rare group (x < 10^-4)

table(msFNA.merged.tcr@meta.data$cloneType)

table(msFNA.merged.tcr@meta.data$ExpansionType)

# check if this frequencies are accurate
ms1Clonotypes = msFNA.merged.tcr@meta.data[msFNA.merged.tcr@meta.data$sampleName=="MS002",]

ctRelFreqs <- sapply(ms1Clonotypes$CTstrict,function(x){
  if(is.na(x)){
    return(NA)
  }else{
    totalNumOfclones = sum(table(ms1Clonotypes$CTstrict)) # not all cells are t cells and/or have tcrs in the data, so we only use total number of tcrs as denominator
    xfreq = sum(ms1Clonotypes$CTstrict==x,na.rm=T)
    clrelfreq = xfreq/totalNumOfclones
    as.numeric(clrelfreq)
  }
  
})

# actual frequencies of clones
ctFreq <- sapply(ms1Clonotypes$CTstrict,function(x){
  if(is.na(x)){
    return(NA)
  }else{
    xfreq = sum(ms1Clonotypes$CTstrict==x,na.rm=T)
    as.numeric(xfreq)
  }
  
})

sum(ms1Clonotypes$Frequency==ctRelFreqs,na.rm=T) 
# so the frequencies calculated by us and scRepertoire are the same
# I also checked that those with highest frequencies were in 10 or more cells
# We will define them as Expanded clones.

# Here the lowest freq we can get in MS002 is 1/5135, as we only  have 5135 clones (which is 0.000194742, and in medium category). So our definition
# is dependent on the number of clones we have. This allows us to divide the clones only into two categories in all samples (expanded and medium)
# Even with this two category definition: we see ebv predicted expanded clones only in MS patients and not controls, also for HCV
# But looking at the actual count frequencies can allow us to divide the clones into four categories (table(ms1Clonotypes$FreqCount))

# We set the predicted antigen species and genes to Unknown if no predication was done by tcrdist3
msFNA.merged.tcr@meta.data$predictedSpecificityVDJdbP[is.na(msFNA.merged.tcr@meta.data$predictedSpecificityVDJdbP)] = "Unknown"
msFNA.merged.tcr@meta.data$predictedEpitopeGeneVDJdbP[is.na(msFNA.merged.tcr@meta.data$predictedEpitopeGeneVDJdbP)] = "Unknown"


## Add expansionType and frequency count meta data to the seurat object
msFNA.merged.tcr@meta.data$ExpansionType = NA
msFNA.merged.tcr@meta.data$FreqCount = NA

for(sname in unique(msFNA.merged.tcr@meta.data$sampleName)){
  snameClonotypes = msFNA.merged.tcr@meta.data[msFNA.merged.tcr@meta.data$sampleName==sname,]
  
  #redefineExpType
  ctRelFreqs <- sapply(snameClonotypes$Frequency,function(x){
    if(is.na(x)){
      return(NA)
    }else{
      if(x > 0.01){
        extype="HyperExpanded"
      }else if(x > 0.001 & x <= 0.01){
        extype="expanded"
      }else if(x > 0.0001 & x <= 0.001){
        extype="medium"
      }else if(x <= 0.0001){
        extype="small"
      }
      return(extype)
    }
    
  })
  
  snameClonotypes$ExpansionType = ctRelFreqs
  
  # Actual frequencies
  ctFreqs <- sapply(snameClonotypes$CTstrict,function(x){
    if(is.na(x)){
      return(NA)
    }else{
      xfreq = sum(snameClonotypes$CTstrict==x,na.rm=T)
      as.numeric(xfreq)
    }
    
  })
  
  snameClonotypes$FreqCount = ctFreqs
  
  # assign the values
  msFNA.merged.tcr@meta.data$ExpansionType[match(rownames(snameClonotypes),rownames(msFNA.merged.tcr@meta.data))] <- snameClonotypes$ExpansionType
  msFNA.merged.tcr@meta.data$FreqCount[match(rownames(snameClonotypes),rownames(msFNA.merged.tcr@meta.data))] <- snameClonotypes$FreqCount
  
}


# Now we can make evaluations and comparisons by expansion.
msFNA.merged.tcr.tcells = msFNA.merged.tcr@meta.data[!is.na(msFNA.merged.tcr@meta.data$barcode),]

specByExpansionType = as.data.frame(table(msFNA.merged.tcr.tcells$ExpansionType,msFNA.merged.tcr.tcells$sampleGroup,msFNA.merged.tcr.tcells$predictedSpecificityVDJdbP))
colnames(specByExpansionType) <- c("Expansion","sampleGroup","Target","Freq")
specByExpansionType2 = specByExpansionType[specByExpansionType$Expansion=="expanded",]

g <- ggplot(data = specByExpansionType2, aes(x = Target, y = Freq,color =sampleGroup))  
g + geom_point() + 
  scale_colour_manual(values=c("lightblue","Brown")) + 
  theme(legend.position = "top") +  theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  ylab("Frequency of expanded clonotypes")

ggsave("MS5_TCR_specificities_and_expandedClones.pdf", height = 210, width = 320, units = "mm")


## % of cells from expanded clonotypes in each sample

# Conclusion: No statistically significant increase in % cells from expanded clones in patients, but overall tendency for increased mean in % of expanded clonotypes 

ExpandedClonesPerSample = as.data.frame(table(msFNA.merged.tcr.tcells$sampleName,msFNA.merged.tcr.tcells$ExpansionType))
colnames(ExpandedClonesPerSample) <- c("sampleName","Expansion","Freq")
ExpandedClonesPerSample$perFreq <- ExpandedClonesPerSample
table(msFNA.merged.tcr.tcells$sampleName)

msFNA.merged.tcr.tcells.expanded = msFNA.merged.tcr.tcells[msFNA.merged.tcr.tcells$ExpansionType=="expanded",]

pExpByExpansionType = c(0,0,0,0,0,0,0,0,0)
names(pExpByExpansionType) <- names(table(msFNA.merged.tcr.tcells$sampleName))

pExpByExpansion = table(msFNA.merged.tcr.tcells.expanded$sampleName)/table(msFNA.merged.tcr.tcells$sampleName)[names(table(msFNA.merged.tcr.tcells.expanded$sampleName))]

pExpByExpansionType[names(pExpByExpansion)] <- pExpByExpansion

pExpByExpansionType <- pExpByExpansionType * 100
pExpByExpansionType <- as.data.frame(pExpByExpansionType)
pExpByExpansionType$sampleGroup = sapply(rownames(pExpByExpansionType),function(jj) ifelse(grepl("CTRL",jj),"CTRL","MS"))

ggplot(pExpByExpansionType, aes(x=sampleGroup, y=pExpByExpansionType,fill=sampleGroup)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:2])+
  geom_label(
    label=rownames(pExpByExpansionType),label.padding = unit(0.08, "lines"), 
    nudge_x = 0.30,nudge_y = -0.05,
    check_overlap = T,size=1
  ) +
  stat_compare_means(aes(group = sampleGroup), label = "p.format", method = "t.test") + 
  theme_classic() + theme(legend.position = "none") +
  ylab("% cells from Expanded clonotypes")

ggsave("MS5_percent_of_expanded_clonotypes_labeled_forPaper.pdf", width = 2, height = 4)

t.test(pExpByExpansionType$pExpByExpansionType[1:3],pExpByExpansionType$pExpByExpansionType[4:9])


## break expanded Tcells by cell type (after review comments)

msFNA.merged.tcr.tcells.expanded = msFNA.merged.tcr.tcells[msFNA.merged.tcr.tcells$ExpansionType=="expanded",]

expandedCellTypesBysample = table(msFNA.merged.tcr.tcells.expanded$rnaclustersvsAdtManAnn2MaxAnn2,msFNA.merged.tcr.tcells.expanded$sampleName)

pheatmap::pheatmap(expandedCellTypesBysample[c("CD4other.2","Effector CD4","Memory CD8"),],scale="none")


msFNA.merged.tcr.tcells.expanded.memCD8s = msFNA.merged.tcr.tcells.expanded[msFNA.merged.tcr.tcells.expanded$rnaclustersvsAdtManAnn2MaxAnn2 %in% c("Memory CD8"),]

msFNA.merged.tcr.tcells.memCD8s <- msFNA.merged.tcr.tcells[msFNA.merged.tcr.tcells$rnaclustersvsAdtManAnn2MaxAnn2 %in% c("Memory CD8"),]

pExpByExpansionType_memCD8s = c(0,0,0,0,0,0,0,0,0)
names(pExpByExpansionType_memCD8s) <- names(table(msFNA.merged.tcr.tcells$sampleName))

pExpByExpansion = table(msFNA.merged.tcr.tcells.expanded.memCD8s$sampleName)/table(msFNA.merged.tcr.tcells.memCD8s$sampleName)[names(table(msFNA.merged.tcr.tcells.expanded.memCD8s$sampleName))]

pExpByExpansionType_memCD8s[names(pExpByExpansion)] <- pExpByExpansion

pExpByExpansionType_memCD8s <- pExpByExpansionType_memCD8s * 100
pExpByExpansionType_memCD8s <- as.data.frame(pExpByExpansionType_memCD8s)
pExpByExpansionType_memCD8s$sampleGroup = sapply(rownames(pExpByExpansionType_memCD8s),function(jj) ifelse(grepl("CTRL",jj),"CTRL","MS"))

ggplot(pExpByExpansionType_memCD8s, aes(x=sampleGroup, y=pExpByExpansionType_memCD8s,fill=sampleGroup)) + 
  geom_boxplot(outlier.colour = NA) +
  geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
  scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:2])+
  geom_label(
    label=rownames(pExpByExpansionType_memCD8s),label.padding = unit(0.08, "lines"), 
    nudge_x = 0.30,nudge_y = -0.05,
    check_overlap = T,size=1
  ) +
  stat_compare_means(aes(group = sampleGroup), label = "p.format", method = "t.test") + 
  theme_classic() + theme(legend.position = "none") +
  ylab("% mem CD8 cells from Expanded mem CD8 clonotypes")

ggsave("MS5_percent_of_expanded_memCD8_clonotypes_labeled_forPaper.pdf", width = 2, height = 4)

t.test(pExpByExpansionType_memCD8s$pExpByExpansionType_memCD8s[1:3],pExpByExpansionType_memCD8s$pExpByExpansionType_memCD8s[4:9])

median(pExpByExpansionType_memCD8s$pExpByExpansionType_memCD8s[4:9])/median(pExpByExpansionType_memCD8s$pExpByExpansionType_memCD8s[1:3])




## % cells from expanded MS, EBV targetting clonotypes in each sample
# Conclusion: no statistical difference in % of expanded clones that target EBV or CMV in MS vs Controls. 

expEBVCMVidx = which(msFNA.merged.tcr.tcells$ExpansionType=="expanded" & !is.na(msFNA.merged.tcr.tcells$predictedSpecificityVDJdbP) & (msFNA.merged.tcr.tcells$predictedSpecificityVDJdbP=="CMV" | msFNA.merged.tcr.tcells$predictedSpecificityVDJdbP=="EBV"))
msFNA.merged.tcr.tcells.expandedEBV_CMV = msFNA.merged.tcr.tcells[expEBVCMVidx,]

pExpByExpansionTypeEBVCMV = c(0,0,0,0,0,0,0,0,0)
names(pExpByExpansionTypeEBVCMV) <- names(table(msFNA.merged.tcr.tcells$sampleName))

pExpByExpansionEBCM = table(msFNA.merged.tcr.tcells.expandedEBV_CMV$sampleName)/table(msFNA.merged.tcr.tcells$sampleName)[names(table(msFNA.merged.tcr.tcells.expandedEBV_CMV$sampleName))]

pExpByExpansionTypeEBVCMV[names(pExpByExpansionEBCM)] <- pExpByExpansionEBCM

pExpByExpansionTypeEBVCMV <- pExpByExpansionTypeEBVCMV * 100
pExpByExpansionTypeEBVCMV <- as.data.frame(pExpByExpansionTypeEBVCMV)
pExpByExpansionTypeEBVCMV$sampleGroup = sapply(rownames(pExpByExpansionTypeEBVCMV),function(jj) ifelse(grepl("CTRL",jj),"CTRL","MS"))

ggplot(pExpByExpansionTypeEBVCMV, aes(x=sampleGroup, y=pExpByExpansionTypeEBVCMV,fill=sampleGroup)) + geom_boxplot() +
  geom_jitter(shape=16, position=position_jitter(0.2)) + 
  scale_colour_manual(values=c("lightblue","Brown")) + 
  theme(legend.position = "top") +  theme_minimal() +
  ylab("% cells from ebv or cmv Expanded clonotypes")

ggsave("/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/TCRanalysis_Results/dbScanResults/MS5_percent_of_EBVCMVexpanded_clonotypes.pdf", height = 210, width = 320, units = "mm")

t.test(pExpByExpansionTypeEBVCMV$pExpByExpansionTypeEBVCMV[1:3],pExpByExpansionTypeEBVCMV$pExpByExpansionTypeEBVCMV[4:9])



## number of expanded clones in sample groups and samples
# What we looked at above is the %cells from expanded clonotypes.  Here we compare the actual number of clonotypes between sample groups. 
# Conclusion: no statistical difference in number of expanded clonotypes between the two
expandedClvsSampleGroups = table(msFNA.merged.tcr.tcells.expanded$CTstrict,msFNA.merged.tcr.tcells.expanded$sampleGroup)
allClonesvsSampleGroups = table(msFNA.merged.tcr.tcells$CTstrict,msFNA.merged.tcr.tcells$sampleGroup)
expandedClonesPer100k = (colSums(expandedClvsSampleGroups > 0) / colSums(allClonesvsSampleGroups > 0)) * 100000

pdf(file="/scratch/project_2005392/JoonaDawitAnalysisResults/combinedAnalysis/TCRanalysis_Results/dbScanResults/MS5_numberOfexpandedClonotypesBySampleGroupPer100k.pdf")
barplot(expandedClonesPer100k,ylim=c(0,60),ylab="#expanded Clonotypes per 100k Clonotypes")
dev.off()

expandedClvsSamples = table(msFNA.merged.tcr.tcells.expanded$CTstrict,msFNA.merged.tcr.tcells.expanded$sampleName)
allClonesvsSamples = table(msFNA.merged.tcr.tcells$CTstrict,msFNA.merged.tcr.tcells$sampleName)
colSums(expandedClvsSamples > 0)
expandedClonesSamplesPer100k = (colSums(expandedClvsSamples > 0) / colSums(allClonesvsSamples > 0)[names(colSums(expandedClvsSamples > 0))]) * 100000
barplot(expandedClonesSamplesPer100k,ylim=c(0,200),ylab="#expanded Clonotypes per 100k Clonotypes")

expandedClonesSamplesPer100kToComp = c(0,0,0,0,0,0,0,0,0)
names(expandedClonesSamplesPer100kToComp) <- names(table(msFNA.merged.tcr@meta.data$sampleName))
expandedClonesSamplesPer100kToComp[names(expandedClonesSamplesPer100k)] <- expandedClonesSamplesPer100k
t.test(expandedClonesSamplesPer100kToComp[4:9],expandedClonesSamplesPer100kToComp[1:3])


expandedClvsCellAnn = table(msFNA.merged.tcr.tcells.expanded$CTstrict,msFNA.merged.tcr.tcells.expanded$rnaclustersvsAdtManAnn2MaxAnn2)
colSums(expandedClvsCellAnn > 0)


# expanded cells by sample
ggplot(msFNA.merged.tcr.tcells.expanded, aes(x=sampleName)) + geom_bar() + 
  geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, colour = "white")
ggsave("MS5_expandedCellsBySamples.pdf", height = 210, width = 320, units = "mm")

# we also removed the NK cells cells with TCR expression
msFNA.merged.tcr.tcells.expanded.tcellsOnly <- msFNA.merged.tcr.tcells.expanded[! msFNA.merged.tcr.tcells.expanded$rnaclustersvsAdtManAnn2MaxAnn2 %in% c("GC B cells","Memory B cells","Monocytes","Myeloid","Naive B cells","Plasmablasts","NK cells"),]

ggplot(msFNA.merged.tcr.tcells.expanded.tcellsOnly, aes(x=sampleName,fill=rnaclustersvsAdtManAnn2MaxAnn2)) + geom_bar() + 
  geom_text(aes(label = ..count..,vjust = 1), position = "stack",stat = "count", colour = "white") + labs(fill = "Annotated cell type") + 
  theme_classic() + theme(legend.position = "right") + ylab("Number of cells expressing expanded TCR clonotypes")
ggsave("MS5_expandedCellsBySamplesByCellType_forPaper.pdf", height = 210, width = 320, units = "mm")

ggplot(msFNA.merged.tcr.tcells.expanded.tcellsOnly, aes(x=sampleName,fill=predictedSpecificityVDJdbP)) + geom_bar(width=.5, position = "fill") + 
  theme_classic() + theme(legend.position = "right") + ylab("Predicted specificity in expanded clonotypes") 
ggsave("MS5_expandedCellsBySamplesBySpecificity2_forPaper.pdf", height = 210, width = 320, units = "mm")


msFNA.merged.tcr.tcells.expanded.tcellsOnly$predictedEpitopeGeneSpeciesVDJdbP <- paste(msFNA.merged.tcr.tcells.expanded.tcellsOnly$predictedEpitopeGeneVDJdbP,
                                                                                       msFNA.merged.tcr.tcells.expanded.tcellsOnly$predictedSpecificityVDJdbP,
                                                                                       sep="-")
msFNA.merged.tcr.tcells.expanded.tcellsOnly$predictedEpitopeGeneSpeciesVDJdbP[msFNA.merged.tcr.tcells.expanded.tcellsOnly$predictedEpitopeGeneSpeciesVDJdbP=="NA-NA"] <- "NA"

ggplot(msFNA.merged.tcr.tcells.expanded.tcellsOnly, aes(x=sampleName,fill=predictedEpitopeGeneSpeciesVDJdbP)) + geom_bar(width=.5, position = "stack") + labs(fill = "Predicted Epitope Genes and species") + 
  theme_classic() + theme(legend.position = "right") + ylab("Predicted specificity in expanded clonotypes") 
ggsave("MS5_expandedCellsBySamplesBySpecificityGenes_forPaper.pdf", height = 210, width = 320, units = "mm")


# split specificity predictions by cell type

ggplot(msFNA.merged.tcr.tcells.expanded.tcellsOnly, aes(x=sampleName,fill=predictedEpitopeGeneSpeciesVDJdbP)) + geom_bar(width=.5, position = "stack") + labs(fill = "Predicted Epitope Genes and species") + 
  theme_classic() + theme(legend.position = "right") + ylab("Predicted specificity in expanded clonotypes") + 
  facet_wrap(vars(rnaclustersvsAdtManAnn2MaxAnn2),ncol=2)

ggsave("MS5_expandedCellTypesBySamplesBySpecificityGenes_forPaper.pdf", height = 210, width = 320,  units = "mm")


ggplot(msFNA.merged.tcr.tcells.expanded.tcellsOnly, aes(x=sampleGroup,fill=predictedEpitopeGeneSpeciesVDJdbP)) + geom_bar(width=.5, position = "stack") + labs(fill = "Predicted Epitope Genes and species") + 
  theme_classic() + theme(legend.position = "right") + ylab("Predicted specificity in expanded clonotypes") + 
  facet_wrap(vars(rnaclustersvsAdtManAnn2MaxAnn2))

ggsave("MS5_expandedCellTypesByGroupsBySpecificityGenes_forPaper.pdf", height = 210, width = 320, units = "mm")


predictedSpecificitiesInExpanded = table(msFNA.merged.tcr.tcells.expanded.tcellsOnly$predictedEpitopeGeneSpeciesVDJdbP,msFNA.merged.tcr.tcells.expanded.tcellsOnly$sampleName)
pheatmap::pheatmap(predictedSpecificitiesInExpanded,scale="none")


# expanded cells by sample Group
ggplot(msFNA.merged.tcr.tcells.expanded, aes(x=sampleGroup)) + geom_bar() + 
  geom_text(aes(label = ..count..), stat = "count", vjust = 1.5, colour = "white")
ggsave("MS5_expandedCellsBySampleGroup.pdf", height = 210, width = 320, units = "mm")

# expanded cells by sample Group by samples
ggplot(msFNA.merged.tcr.tcells.expanded, aes(x=sampleGroup,fill=sampleName)) + geom_bar() 
ggsave("MS5_expandedCellsBySampleGroupBySamples.pdf", height = 210, width = 320, units = "mm")




## Evaluate TCR sharing between t cell subsets, and between MS and CTRLs ####

msFNA.merged.tcr.tcells = msFNA.merged.tcr@meta.data[!is.na(msFNA.merged.tcr@meta.data$barcode),]

tcrsAnnotatedAsTcells = c("Effector CD4", "Tfh.1", "Treg ", "Naive CD4", "Tfr", "Intermediate CD4", "Tfh.2", "Naive CD8", "Memory CD8",
                          "CD4other.2", "CD4other.1")

msFNA.merged.tcr.tcellsOnly <- msFNA.merged.tcr.tcells[msFNA.merged.tcr.tcells$rnaclustersvsAdtManAnn2MaxAnn2 %in% tcrsAnnotatedAsTcells, ]

TCRoverlapBetweenTcellSubsets_num = NULL
TCRoverlapBetweenTcellSubsets_overlapRate = NULL
TCRoverlapBetweenTcellSubsets_ClonalityOfSharedClones = NULL # as percent of non-unique clones or copies
TCRoverlapBetweenTcellSubsets_Clonality = NULL # as percent of non-unique clones or copies


sams = unique(msFNA.merged.tcr.tcellsOnly$sampleName)

follicularTcells = c("Tfh.1","Tfh.2","Tfr")

for(sam in sams){
  msFNA.merged.tcr.tcellsOnly.sam = msFNA.merged.tcr.tcellsOnly[msFNA.merged.tcr.tcellsOnly$sampleName==sam,]
  
  fTcellByOtherTcell_Overlaps = c()
  fTcellByOtherTcell_Nums = c()
  clonality_shared = c()
  clonalitySubset = c()
  
  
  for(fTcells in follicularTcells){
    
    others_tcrsAnnotatedAsTcells =  tcrsAnnotatedAsTcells
    
    fTcells_tcrs_sam = msFNA.merged.tcr.tcellsOnly.sam$CTstrict[msFNA.merged.tcr.tcellsOnly.sam$rnaclustersvsAdtManAnn2MaxAnn2==fTcells]
    
    # clonality in fTcell subset
    clonalityInSamfTcells = (length(fTcells_tcrs_sam) - length(unique(fTcells_tcrs_sam))) / length(fTcells_tcrs_sam) # clonality 
    names(clonalityInSamfTcells) <- fTcells
    clonalitySubset = c(clonalitySubset,clonalityInSamfTcells * 100)
    
    for(otherCells in others_tcrsAnnotatedAsTcells){
      otherCells_tcrs_sam = msFNA.merged.tcr.tcellsOnly.sam$CTstrict[msFNA.merged.tcr.tcellsOnly.sam$rnaclustersvsAdtManAnn2MaxAnn2==otherCells]
      
      # overlap calculation
      numberOfClonesSam1 = length(unique(fTcells_tcrs_sam))
      numberOfClonesSam2 = length(unique(otherCells_tcrs_sam))
      numberOfSharedClones = length(intersect(fTcells_tcrs_sam,otherCells_tcrs_sam))
      
      ov = 2 * (numberOfSharedClones / (numberOfClonesSam1 + numberOfClonesSam2))
      
      names(ov) <- paste(fTcells,otherCells,sep="_")
      names(numberOfSharedClones) <- paste(fTcells,otherCells,sep="_")
      
      fTcellByOtherTcell_Overlaps <- c(fTcellByOtherTcell_Overlaps,ov)
      fTcellByOtherTcell_Nums <- c(fTcellByOtherTcell_Nums,numberOfSharedClones)
      
      # clonalityWithinEachSubset of shared clones
      if(numberOfSharedClones==0){
        clonalityOfShared = 0
      }else{
        
        fTcells_tcrs_sam_shared = fTcells_tcrs_sam[fTcells_tcrs_sam %in% intersect(fTcells_tcrs_sam,otherCells_tcrs_sam)]
        otherCells_tcrs_sam_shared = otherCells_tcrs_sam[otherCells_tcrs_sam %in% intersect(fTcells_tcrs_sam,otherCells_tcrs_sam)]
        
        clonalitySam1 = (length(fTcells_tcrs_sam_shared) - length(unique(fTcells_tcrs_sam_shared))) / length(fTcells_tcrs_sam_shared) # clonality 
        clonalitySam2 = (length(otherCells_tcrs_sam_shared) - length(unique(otherCells_tcrs_sam_shared))) / length(otherCells_tcrs_sam_shared) # clonality 
        
        allShared = c(fTcells_tcrs_sam_shared,otherCells_tcrs_sam_shared)
        #allShared = c(fTcells_tcrs_sam,otherCells_tcrs_sam) # how clonal are the two subsets together.
        
        #clonalityOfShared = (length(allShared) - length(unique(allShared))) / max(length(fTcells_tcrs_sam),length(otherCells_tcrs_sam))
        
        clonalityOfShared = mean(c(clonalitySam1,clonalitySam2))
      }
      
      names(clonalityOfShared) <- paste(fTcells,otherCells,sep="_")
      clonality_shared = c(clonality_shared,clonalityOfShared * 100)
      
    }
    
  }
  
  TCRoverlapBetweenTcellSubsets_num <- cbind(TCRoverlapBetweenTcellSubsets_num,fTcellByOtherTcell_Nums)
  TCRoverlapBetweenTcellSubsets_overlapRate <- cbind(TCRoverlapBetweenTcellSubsets_overlapRate,fTcellByOtherTcell_Overlaps)
  TCRoverlapBetweenTcellSubsets_ClonalityOfSharedClones <- cbind(TCRoverlapBetweenTcellSubsets_ClonalityOfSharedClones,clonality_shared)
  TCRoverlapBetweenTcellSubsets_Clonality <- cbind(TCRoverlapBetweenTcellSubsets_Clonality,clonalitySubset)
  
}


colnames(TCRoverlapBetweenTcellSubsets_num) <- sams
colnames(TCRoverlapBetweenTcellSubsets_overlapRate) <- sams
colnames(TCRoverlapBetweenTcellSubsets_ClonalityOfSharedClones) <- sams
colnames(TCRoverlapBetweenTcellSubsets_Clonality) <- sams

# TCRsubset overlapping t test between groups


TCRsubsetSharingPvals = apply(TCRoverlapBetweenTcellSubsets_overlapRate,1,function(x){
  g1 = x[grepl("MS",names(x))]
  g2 = x[grepl("CTRL",names(x))]
  
  toR = tryCatch(expr = {return(t.test(g1,g2)$p.value)},error=function(e){return(1)})
  toR
})

TCRsubsetClonalitySharedPvals = apply(TCRoverlapBetweenTcellSubsets_ClonalityOfSharedClones,1,function(x){
  g1 = x[grepl("MS",names(x))]
  g2 = x[grepl("CTRL",names(x))]
  
  toR = tryCatch(expr = {return(t.test(g1,g2)$p.value)},error=function(e){return(1)})
  toR
})

TCRsubsetClonalityPvals = apply(TCRoverlapBetweenTcellSubsets_Clonality,1,function(x){
  g1 = x[grepl("MS",names(x))]
  g2 = x[grepl("CTRL",names(x))]
  
  toR = tryCatch(expr = {return(t.test(g1,g2)$p.value)},error=function(e){return(1)})
  toR
})

#
TCRoverlapBetweenTcellSubsets_overlapRate_pval <- cbind(TCRoverlapBetweenTcellSubsets_overlapRate,TCRsubsetSharingPvals)
TCRoverlapBetweenTcellSubsets_ClonalityOfSharedClones_pval <- cbind(TCRoverlapBetweenTcellSubsets_ClonalityOfSharedClones,TCRsubsetClonalitySharedPvals)
TCRoverlapBetweenTcellSubsets_Clonality_pval <- cbind(TCRoverlapBetweenTcellSubsets_Clonality,TCRsubsetClonalityPvals)


write.csv(TCRoverlapBetweenTcellSubsets_num,"TCRoverlapBetweenTcellSubsets_num.csv")
write.csv(TCRoverlapBetweenTcellSubsets_overlapRate_pval,"TCRoverlapBetweenTcellSubsets_overlapRate_pval.csv")
write.csv(TCRoverlapBetweenTcellSubsets_ClonalityOfSharedClones_pval,"TCRoverlapBetweenTcellSubsets_ClonalityOfSharedClones_pval.csv")
write.csv(TCRoverlapBetweenTcellSubsets_Clonality_pval,"TCRoverlapBetweenTcellSubsets_Clonality_pval.csv")


## heatmaps of TCR subset sharing, for follicular T cells

# TCRsubset overlapRate, first remove same pairs that bias the plots, and those that are the repeats of the same pairs
removePairsAndRepeats = c("Tfh.1_Tfh.1","Tfh.2_Tfh.2","Tfr_Tfr","Tfh.2_Tfh.1","Tfr_Tfh.1","Tfr_Tfh.2")
TCRoverlapBetweenTcellSubsets_overlapRate_toPlot = TCRoverlapBetweenTcellSubsets_overlapRate[!rownames(TCRoverlapBetweenTcellSubsets_overlapRate) %in% removePairsAndRepeats,]

pdf(file="TCRoverlapBetweenTcellSubsets_overlapRate_toPlot_heatmap.pdf")
pheatmap::pheatmap(TCRoverlapBetweenTcellSubsets_overlapRate_toPlot)
dev.off()

# TCRsubset overlap counts
TCRoverlapBetweenTcellSubsets_num_toPlot = TCRoverlapBetweenTcellSubsets_num[!rownames(TCRoverlapBetweenTcellSubsets_num) %in% removePairsAndRepeats,]

pdf(file="TCRoverlapBetweenTcellSubsets_num_heatmap.pdf")
pheatmap::pheatmap(TCRoverlapBetweenTcellSubsets_num_toPlot)
dev.off()



## Plot results of detected TCR types ####

library(reshape2)
library(ggsignif)
# dot plots of TCR specificities in each sample

longnHitsPerPathTypeFull <- melt(data = nHitsPerPathTypeFull[-1,1:9] + 1,variable.name = "samples", 
                                 value.name = "nSpecificTCRs")

longnHitsPerPathTypeFull$sampleGroup <- sapply(longnHitsPerPathTypeFull$Var2,function(x) ifelse(grepl("CTRL",x),"CTRL","MS"))
longnHitsPerPathTypeFull$pathology <- longnHitsPerPathTypeFull$Var1
longnHitsPerPathTypeFull$nSpecificTCRsLog <- log10(longnHitsPerPathTypeFull$nSpecificTCRs)

g <- ggplot(data = longnHitsPerPathTypeFull, aes(x = pathology, y = nSpecificTCRsLog,fill =sampleGroup))  
g + geom_boxplot(aes(colour = sampleGroup),outlier.shape = NA) + 
  geom_point(aes(color= sampleGroup, shape=sampleGroup), alpha = 0.9, size=1.5,position = position_jitterdodge(jitter.width = 0.1)) + 
  scale_fill_manual(values=c("lightblue","Brown")) +
  #scale_colour_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:2]) + 
  theme_classic() + theme(legend.position = "none", axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("# specific TCRs per 1000 T-cell clonotypes")
# theme(legend.position = "top") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +


ggsave("dbScanResults/VDJdb_MS5_TCR_specificities_profile.pdf", height = 210, width = 320, units = "mm")


## plot the TCR specificity profile for all T cell subsets

for(cType in names(nHitsPerPathTypeList)){
  
  longnHitsPerPathTypeFull <- melt(data = nHitsPerPathTypeList[[cType]][-1,1:9] + 1,variable.name = "samples", 
                                   value.name = "nSpecificTCRs")
  
  longnHitsPerPathTypeFull$sampleGroup <- sapply(longnHitsPerPathTypeFull$Var2,function(x) ifelse(grepl("CTRL",x),"CTRL","MS"))
  longnHitsPerPathTypeFull$pathology <- longnHitsPerPathTypeFull$Var1
  longnHitsPerPathTypeFull$nSpecificTCRsLog <- log10(longnHitsPerPathTypeFull$nSpecificTCRs)
  
  g <- ggplot(data = longnHitsPerPathTypeFull, aes(x = pathology, y = nSpecificTCRsLog,color =sampleGroup))  
  g + geom_boxplot(outlier.shape=NA,position = position_dodge(width = 0.6),alpha = 0.8) + 
    geom_point(aes(color= sampleGroup, shape=sampleGroup), alpha = 0.7, size=1.5,position = position_jitterdodge(jitter.width = 0.1)) + 
    #geom_dotplot(binaxis='y', stackdir = "center", position = "dodge",dotsize=0.5,alpha=0.4) +
    scale_color_manual(values=c("lightblue","Brown")) + 
    #scale_colour_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:2]) + 
    theme_minimal() + 
    theme(legend.position = "top", axis.text.x = element_text(angle = 45, hjust = 1)) + 
    ylab("# specific TCRs per 1000 T-cells (log10)")
  # theme(legend.position = "top") + theme_minimal() + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  
  
  ggsave(paste0("dbScanResults/",cType,"_VDJdb_MS5_paired_TCR_specificities_profile_finalAnnotation_afterHLAFiltering.pdf"), height = 210, width = 320, units = "mm")
  
  
}

# print results for T cell subpopulations that show significantly different TCR specificity between the groups
for(ctype in names(nHitsPerPathTypeList)){
  
  print(ctype)
  pvalCol = which(grepl("value",colnames(nHitsPerPathTypeList[[ctype]])))
  
  sigVals = which(nHitsPerPathTypeList[[ctype]][,pvalCol] < 0.09)
  
  if(length(sigVals) > 0){
    print(nHitsPerPathTypeList[[ctype]][sigVals,,drop=F])
  }
  
}



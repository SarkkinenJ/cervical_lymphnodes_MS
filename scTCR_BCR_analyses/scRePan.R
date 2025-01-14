source("TCR_functions.R")

# for now reps should be a list provided by scRepertoire's combineTCR function. It is a named list of samples with each element holding the repertoire data per sample
# e.g combined_MS5_TCRcontigList comes from MSFNA_MergedData_TCRanalysis.R after reading in all scTCR data for each sample

reps = combined_MS5_TCRcontigList


#............Inter-individual overlap analysis.............#####

# estimate inter-individual overlap using O= 2 x [c /(a+b) ] 

MSFNA_msTCR_InterIndivualOverlapNTpaired <- estimateInterIndividualOverlaps(reps[grepl("MS",names(reps))],downsampleSize=downSampleTo,useSampling=T,useSeq="CTstrict",nResamples=10)
MSFNA_ctrlTCR_InterIndivualOverlapNTpaired <- estimateInterIndividualOverlaps(reps[grepl("CTRL",names(reps))],downsampleSize=downSampleTo,useSampling=T,useSeq="CTstrict",nResamples=10)
t.test(as.numeric(MSFNA_msTCR_InterIndivualOverlapNTpaired[,2]),as.numeric(MSFNA_ctrlTCR_InterIndivualOverlapNTpaired[,2])) # a tendency towards increased sharing among MS patients


MSFNA_msTCR_InterIndivualOverlapAApaired <- estimateInterIndividualOverlaps(reps[grepl("MS",names(reps))],downsampleSize=downSampleTo,useSampling=T,useSeq="CTaa",nResamples=10)
MSFNA_ctrlTCR_InterIndivualOverlapAApaired <- estimateInterIndividualOverlaps(reps[grepl("CTRL",names(reps))],downsampleSize=downSampleTo,useSampling=T,useSeq="CTaa",nResamples=10)

t.test(as.numeric(MSFNA_msTCR_InterIndivualOverlapAApaired[,2]),as.numeric(MSFNA_ctrlTCR_InterIndivualOverlapAApaired[,2])) # a tendency towards increased sharing among MS patients


MSFNA_msTCR_InterIndivualOverlapAAalpha <- estimateInterIndividualOverlaps(reps[grepl("MS",names(reps))],downsampleSize=downSampleTo,useSampling=T,useSeq="cdr3_aa1",nResamples=10)
MSFNA_ctrlTCR_InterIndivualOverlapAAalpha <- estimateInterIndividualOverlaps(reps[grepl("CTRL",names(reps))],downsampleSize=downSampleTo,useSampling=T,useSeq="cdr3_aa1",nResamples=10)

t.test(as.numeric(MSFNA_msTCR_InterIndivualOverlapAAalpha[,2]),as.numeric(MSFNA_ctrlTCR_InterIndivualOverlapAAalpha[,2])) # a tendency towards increased sharing among MS patients


MSFNA_msTCR_InterIndivualOverlapAAbeta <- estimateInterIndividualOverlaps(reps[grepl("MS",names(reps))],downsampleSize=downSampleTo,useSampling=T,useSeq="cdr3_aa2",nResamples=10)
MSFNA_ctrlTCR_InterIndivualOverlapAAbeta <- estimateInterIndividualOverlaps(reps[grepl("CTRL",names(reps))],downsampleSize=downSampleTo,useSampling=T,useSeq="cdr3_aa2",nResamples=10)

t.test(as.numeric(MSFNA_msTCR_InterIndivualOverlapAAbeta[,2]),as.numeric(MSFNA_ctrlTCR_InterIndivualOverlapAAbeta[,2])) # a tendency towards increased sharing among MS patients


## plots inter-individial overlap

    # plot diversity dot plot
    interIndOverlapAApaired = data.frame(overlap=c(MSFNA_msTCR_InterIndivualOverlapAApaired[,2],MSFNA_ctrlTCR_InterIndivualOverlapAApaired[,2]),
                                   sampleGroup=c(rep("MS",15),rep("CTRL",3)))
    interIndOverlapAApaired$overlapType = "ab"
    
    interIndOverlapAAalpha = data.frame(overlap=c(MSFNA_msTCR_InterIndivualOverlapAAalpha[,2],MSFNA_ctrlTCR_InterIndivualOverlapAAalpha[,2]),
                                   sampleGroup=c(rep("MS",15),rep("CTRL",3)))
    interIndOverlapAAalpha$overlapType = "a"
    
    
    interIndOverlapAABeta = data.frame(overlap=c(MSFNA_msTCR_InterIndivualOverlapAAbeta[,2],MSFNA_ctrlTCR_InterIndivualOverlapAAbeta[,2]),
                                        sampleGroup=c(rep("MS",15),rep("CTRL",3)))
    interIndOverlapAABeta$overlapType = "b"
    
    interIndOverlapAA = rbind(interIndOverlapAApaired,interIndOverlapAAalpha,interIndOverlapAABeta)
    
    interIndOverlapAA$sampleGroup <- as.factor(interIndOverlapAA$sampleGroup)
    interIndOverlapAA$overlap <- as.numeric(interIndOverlapAA$overlap)
    
    interIndOverlapAA_statp  = compare_means(overlap ~ sampleGroup, interIndOverlapAA,group.by="overlapType", method="t.test",p.adjust.method="BH")
    
    stat.test <- interIndOverlapAA_statp %>%
      mutate(y.position = rep(0.025,3))
    
  for(overlapT in unique(interIndOverlapAA$overlapType)){
    
    interIndOverlapAAD = interIndOverlapAA[interIndOverlapAA$overlapType==overlapT,]
    interIndOverlapAA_statpD = interIndOverlapAA_statp[interIndOverlapAA_statp$overlapType==overlapT,]
    
    # for overlap type ab, all shared clones had one NA chain, so really no sharing so it is set to zero.
    if(overlapT=="ab"){
      interIndOverlapAAD$overlap <- 0
      }
    
    stat.testD <- interIndOverlapAA_statpD %>%
      mutate(y.position = rep(max(interIndOverlapAAD$overlap) + 0.0001,1))
    
    
    g <- ggplot(data = interIndOverlapAAD, aes(x = sampleGroup,y = overlap,fill=sampleGroup))  
    g +  geom_boxplot(outlier.colour = NA) +  
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") + 
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Inter-individual repertoire overlap") + ggtitle(paste0("TCR chain: ",overlapT)) + 
      
      # significance indicated with stars
      #stat_pvalue_manual(stat.test, x = "timePoint", hide.ns = TRUE,label = "p.signif",label.size=6,inherit.aes = FALSE,tip.length = 0)
      
      # significance with pvalues
      stat_pvalue_manual(stat.testD, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)
    
    ggsave(paste0("MS5_scRepAn_interIndividualAAClonotypeOverlap_",overlapT,".pdf"), width = 2, height = 4)
    
  }


#............Diversity analysis after downsampling .............#####

downSampleTo=findSmallestNumberOfReads(reps) # smallest repertoire MS002 and has size of 5258 clonotypes

MSFNA_msTCR_Diversity_downsampled <- estimateDiversity(reps[grepl("MS",names(reps))],downsampleSize=downSampleTo,useSeq="nt",nResamples=100)

MSFNA_ctrlTCR_Diversity_downsampled <- estimateDiversity(reps[grepl("CTRL",names(reps))],downsampleSize=downSampleTo,useSeq="nt",nResamples=100)


t.test(MSFNA_msTCR_Diversity_downsampled,MSFNA_ctrlTCR_Diversity_downsampled) # no difference in diversity

# plot diversity dot plot
diversityValues = data.frame(diversity=c(MSFNA_msTCR_Diversity_downsampled,MSFNA_ctrlTCR_Diversity_downsampled),
                             sampleGroup=c(rep("MS",length(MSFNA_msTCR_Diversity_downsampled)),
                                           rep("CTRL",length(MSFNA_ctrlTCR_Diversity_downsampled)))
                             )

diversityValues$sampleGroup <- as.factor(diversityValues$sampleGroup)

p <- ggplot(diversityValues, aes(x=sampleGroup, y=diversity)) + 
  geom_violin(trim = FALSE)+
  geom_dotplot(binaxis='y', stackdir='center') + labs(y = "Normalized Shannon Entropy") + theme_bw()

p + stat_summary(fun.data="mean_sdl", fun.args = list(mult=1), 
                   geom="crossbar", width=0.5)

ggsave("MS5_scRepAn_diversityAnalysisNtClonotype.pdf", height = 210, width = 320, units = "mm")


## amino acid level diversity
MSFNA_msTCR_Diversity_downsampledAA <- estimateDiversity(reps[grepl("MS",names(reps))],downsampleSize=downSampleTo,useSeq="AA",nResamples=100)

MSFNA_ctrlTCR_Diversity_downsampledAA <- estimateDiversity(reps[grepl("CTRL",names(reps))],downsampleSize=downSampleTo,useSeq="AA",nResamples=100)

t.test(MSFNA_msTCR_Diversity_downsampledAA,MSFNA_ctrlTCR_Diversity_downsampledAA) # no difference in diversity again in AA

# plot diversity dot plot
diversityValuesAA = data.frame(diversity=c(MSFNA_msTCR_Diversity_downsampledAA,MSFNA_ctrlTCR_Diversity_downsampledAA),
                             sampleGroup=c(rep("MS",length(MSFNA_msTCR_Diversity_downsampledAA)),
                                           rep("CTRL",length(MSFNA_ctrlTCR_Diversity_downsampledAA)))
)

diversityValuesAA$sampleGroup <- as.factor(diversityValuesAA$sampleGroup)

diversityValuesAA_statp  = compare_means(diversity ~ sampleGroup, diversityValuesAA,method="t.test",p.adjust.method="BH")

stat.test <- diversityValuesAA_statp %>%
  mutate(y.position = rep(0.99,1))


    g <- ggplot(data = diversityValuesAA, aes(x = sampleGroup,y = diversity,fill=sampleGroup))  
    g + geom_boxplot(outlier.colour = NA) + 
      #scale_fill_viridis(discrete = TRUE, alpha=0.6, option="A") +
      #scale_color_viridis(discrete = TRUE, alpha=0.6, option="A") + 
      #scale_color_manual(values=c("lightblue","Brown")) + 
      scale_fill_manual(values = RColorBrewer::brewer.pal(3, "Accent")[1:3]) + 
      geom_point(position=position_jitterdodge(jitter.width=0), pch=21) +
      theme_classic() + theme(legend.position = "none") + 
      theme(axis.line = element_line(color = "grey70"),plot.title = element_text(hjust = 0.5)) + 
      ylab("Normalized Shannon Entropy") +
      
      #stat_pvalue_manual(stat.test, x = "timePoint", hide.ns = TRUE,label = "p.signif",label.size=6,inherit.aes = FALSE,tip.length = 0)
      
      # significance with pvalues
      stat_pvalue_manual(stat.test, hide.ns = F,label = "p = {scales::pvalue(p)}",inherit.aes = FALSE,tip.length = 0)

    ggsave("MS5_scRepAn_diversityAnalysisAAClonotype.pdf",width = 2, height = 4)
    

ggsave("MS5_scRepAn_diversityAnalysisAAClonotype.pdf", height = 210, width = 320, units = "mm")








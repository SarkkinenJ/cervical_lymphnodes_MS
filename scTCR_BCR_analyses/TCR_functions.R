
findSmallestNumberOfReads <- function(reps){
  
  sampleNames = names(reps)
  sampleSizes = c()
  for(r in names(reps)){
    
    sam1 <- reps[[r]]

    sampleSizes = c(sampleSizes,nrow(sam1)) # this should probably be sum(table(sam1$CTstrict)) as some clonotypes are seen repeatedly
  }
  
  #print(sampleSizes)
  #print(sampleNames[which.min(sampleSizes)])
  
  return(min(sampleSizes))
  
}

estimateOverlap <- function(sam1,sam2,downsampleSize=NULL,useSeq="CTaa",nResamples=10){
  
  
  if(is.null(downsampleSize)){
    resampleSize = min(nrow(sam1),nrow(sam2))
  }else{
    resampleSize = downsampleSize
  }
  
  overlaps = c()
  overlappingSeqs = c()
  
  
  for(i in 1:nResamples){
    
    #sam1 resample

    resampledSam1_idx = sample(1:nrow(sam1),resampleSize,replace=T)
    resampledSam1 = sam1[resampledSam1_idx,]
    
    #sam2 resample

    resampledSam2_idx = sample(1:nrow(sam2),resampleSize,replace=T)
    resampledSam2 = sam2[resampledSam2_idx,]
    
   
    # make selection of clonotypes by paired Nt, AA or single chain Nt or AA
    numberOfClonesSam1 = length(unique(resampledSam1[,useSeq][!is.na(resampledSam1[,useSeq])]))
    numberOfClonesSam2 = length(unique(resampledSam2[,useSeq][!is.na(resampledSam2[,useSeq])]))
    
    numberOfSharedClones = length(intersect(resampledSam1[,useSeq][!is.na(resampledSam1[,useSeq])],
                                            resampledSam2[,useSeq][!is.na(resampledSam2[,useSeq])]))
    
  
    print(numberOfClonesSam1)
    print(numberOfClonesSam2)
    print(numberOfSharedClones)
  
    ov = 2 * (numberOfSharedClones / (numberOfClonesSam1 + numberOfClonesSam2))
    
    overlaps = c(overlaps,ov)
    overlappingSeqs <- c(overlappingSeqs,intersect(resampledSam1[,useSeq],resampledSam2[,useSeq]))
    
  }
  
  #print(overlaps)
  medianOverlap = median(overlaps)
  overlappingSeqs = unique(overlappingSeqs)
  print(overlappingSeqs)
  
  return(list(medianOverlap,overlappingSeqs))
  
}

estimateOverlapNoDownsampling <- function(sam1,sam2,useSeq="CTaa"){
  
  # calculate overlap
  
  numberOfClonesSam1 = length(unique(sam1[,useSeq]))
  numberOfClonesSam2 = length(unique(sam2[,useSeq]))
  numberOfSharedClones = length(intersect(sam1[,useSeq],sam2[,useSeq]))
  
  ov = 2 * (numberOfSharedClones / (numberOfClonesSam1 + numberOfClonesSam2))
  
  
  #print(overlaps)
  
  return(ov)
  
}


estimateInterIndividualOverlaps <- function(reps,useSampling=T,downsampleSize=NULL,useSeq="CTaa",nResamples=10){
  
  temp <- NULL
  nm <- NULL
  overlapSeqs = NULL
  
  numberOfSamples = length(reps)
  
  for(i in 1:(numberOfSamples-1)){
    
    nextIndex = i+1
    
    for(j in nextIndex:numberOfSamples){
      
      print(names(reps)[i])
      print(names(reps)[j])
      
      if(useSampling==T){
        op <- estimateOverlap(reps[[i]],reps[[j]],downsampleSize,useSeq,nResamples)
      }else{
        op <- estimateOverlapNoDownsampling(reps[[i]],reps[[j]],useSeq)
        
      }
      #m <- paste(rep,"vs",r,sep="")
      #nm <- c(nm,m)
      temp <- rbind(temp,op[[1]])
      nm <- rbind(nm,paste(names(reps)[i],names(reps)[j],sep=" Vs "))
      overlapSeqs <- rbind(overlapSeqs,paste0(op[[2]],collapse=":"))
      #overlapSeqs[[nm]] <- op[[2]]
    }
    
  }
  
  toRet = cbind(nm,temp,overlapSeqs)  
  return(toRet)
}


# shannon index: parameter is clonotype frequency
shannon_index=function(f){
  
  tot_sum=sum(f)
  fi=f/tot_sum
  
  first <- fi * log(fi)
  sh <- abs(sum(first))
  #sh_effective=exp(sh) # effective size of repertoire
  sh <- sh/log(length(f)) # normalizing entropy to number of unique TCRs
  sh <- 1-sh # clonality
  return(sh)
}

estimateDiversity <- function(reps,downsampleSize=NULL,useSeq="AA",minCTSize=1,nResamples=10){
  
  entropies=c()
  
  for(rep in reps){
    
    dsTemp <- c()
    
    for(i in 1:nResamples){
      #sam1 resample

      resampledSam1_idx = sample(1:nrow(rep),downsampleSize,replace=T)
      resampledSam1 = rep[resampledSam1_idx,]
        
        if(useSeq=="AA"){
          resampledSam1Count = sort(table(resampledSam1$CTaa),decreasing=T) 
        }else{
          resampledSam1Count = sort(table(resampledSam1$CTstrict),decreasing=T) # nt clones must include same genes nt+genes
        }
        
     
      counts=resampledSam1Count[resampledSam1Count >= minCTSize]
      p_ds <- 1-shannon_index(counts) # get the diversity
      dsTemp <- c(dsTemp,p_ds) 
      
    }
    
    entropies <- c(entropies,median(dsTemp))
    
  }
  
  
  return(entropies)
}

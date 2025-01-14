library(dplyr)
library(Seurat)
library(ggplot2)

library(scRepertoire)
library("reticulate")

# read in scTCR data for samples ####

combined_MS5_TCRcontigList <- readRDS(file = "data/combined_MS5_TCRcontigList.rds")



# run the db scan ####


pd <- import("pandas")
tcrdist <- import("tcrdist")


# VDJdb latest version (30 March, 2022)
# download vdjdb data here and read it in

VDJdblatest = read.csv2("vdjdb_full.txt",sep="\t")
VDJdblatest[VDJdblatest==""] <- NA
VDJdblatest = VDJdblatest[VDJdblatest$species=="HomoSapiens",]


VDJdblatest$did = rownames(VDJdblatest)
VDJdblatest$count = 1
mainCols = c("did","meta.subject.id","antigen.epitope","count","v.alpha","j.alpha","cdr3.alpha","v.beta","j.beta","cdr3.beta")
otherCols = colnames(VDJdblatest)[!colnames(VDJdblatest) %in% mainCols]
newCols = c(mainCols,otherCols)
VDJdblatest <- VDJdblatest[,newCols]
colnames(VDJdblatest)[1:10] <- colnames(dfexample)


#  first column should be id
colnames(VDJdblatest)[1] <- "did"



# records with paired TCR only
VDJdblatest_paired = VDJdblatest[!(is.na(VDJdblatest$v_a_gene) | is.na(VDJdblatest$j_a_gene) | is.na(VDJdblatest$v_b_gene) | is.na(VDJdblatest$j_b_gene)),] 

# TCRb records only
VDJdblatest_tcrb = VDJdblatest[!(is.na(VDJdblatest$v_b_gene) | is.na(VDJdblatest$j_b_gene)),] 


## run the scanning 
numberOfEBVhits <- NULL
numberOfTCRs <- NULL
proportionOfGenesDetectedIndatabase = NULL

missedNums <- NULL


toWrite = "dbScanResults/"

# defailts to paired analysis
args <- commandArgs(trailingOnly = TRUE)
dbtype <- args[1]

if(dbtype=="paired" | is.null(dbtype)){
  dbtype="paired"
  VDJdblatest_db = VDJdblatest_paired
}else{
  dbtype="beta"
  VDJdblatest_db = VDJdblatest_tcrb
}

# the sample to analyze
sampleidArg <- as.numeric(args[2])
combined_MS5_TCRcontigList <- combined_MS5_TCRcontigList[sampleidArg]

print(names(combined_MS5_TCRcontigList))


for(sam in names(combined_MS5_TCRcontigList)){
  sampTCRs <- combined_MS5_TCRcontigList[[sam]]
  cloneCount <- table(sampTCRs$CTstrict)
  sampTCRdata = data.frame(did=rownames(sampTCRs),subject=sampTCRs$sample)
  sampTCRdata$epitope <- "None"
  sampTCRdata$count <- as.integer(cloneCount[sampTCRs$CTstrict])
  sampTCRdata$v_a_gene <- sapply(sampTCRs$TCR1,function(x) ifelse(is.na(unlist(strsplit(x,"\\."))[1]),"None",paste0(unlist(strsplit(x,"\\."))[1],"*01")))
  sampTCRdata$j_a_gene <- sapply(sampTCRs$TCR1,function(x) ifelse(is.na(unlist(strsplit(x,"\\."))[2]),"None",paste0(unlist(strsplit(x,"\\."))[2],"*01")))
  sampTCRdata$cdr3_a_aa <- sapply(sampTCRs$cdr3_aa1,function(x) unlist(strsplit(x,";"))[1])
  sampTCRdata$v_b_gene <- sapply(sampTCRs$TCR2,function(x) ifelse(is.na(unlist(strsplit(x,"\\."))[1]),"None",paste0(unlist(strsplit(x,"\\."))[1],"*01")))
  sampTCRdata$j_b_gene <-  sapply(sampTCRs$TCR2,function(x) ifelse(is.na(unlist(strsplit(x,"\\."))[2]),"None",paste0(unlist(strsplit(x,"\\."))[2],"*01")))
  sampTCRdata$cdr3_b_aa <- sapply(sampTCRs$cdr3_aa2,function(x) unlist(strsplit(x,";"))[1])

  # remove clones that don't have both chains. This has to be revisited
  sampTCRdata <- sampTCRdata[apply(sampTCRdata[,5:10],1, function(x) sum("None" %in% x) == 0 & sum(is.na(x)) == 0),]
  
  numberOfTCRs <- c(numberOfTCRs,nrow(sampTCRdata))
 
  sampTCRdata_andDatabase <- rbind(sampTCRdata,VDJdblatest_db[,colnames(dfexample)])
  sampTCRdata_andDatabase$did <- as.integer(sampTCRdata_andDatabase$did)
  
  # divide the database into chunks of 5000 records and run the tcrdist analysis
  tr <- list()
  
  chunkSize = floor(10000 - nrow(sampTCRdata))
  
  nChunks = ceiling(nrow(VDJdblatest_db)/chunkSize)
  startPt = 1
  
  for(i in 1:nChunks){
    
    
    endr = (startPt + chunkSize)-1
    if(endr > nrow(VDJdblatest_db)){
      endr = nrow(VDJdblatest_db)
    }
    
    print(i)
    print(startPt)
    print(endr)
    print("---")
    
    
    dfChunk = VDJdblatest_db[startPt:endr,]
    startPt = startPt + chunkSize
    
    sampTCRdata_andDatabaseT <- rbind(sampTCRdata,dfChunk[,colnames(dfexample)])
    sampTCRdata_andDatabaseT$id <- 1:nrow(sampTCRdata_andDatabaseT)
     
    
    trTemp = tcrdist$repertoire$TCRrep(cell_df = sampTCRdata_andDatabaseT,organism = 'human', chains = c('alpha', 'beta'),
                                       db_file = 'alphabeta_gammadelta_db.tsv')
    
    tr[[i]] <- trTemp
    
    
  }
  
  res2 = NULL
  nEBVForClone = c()
  
  for(nch in 1:length(tr)){
    
    
    
    samcloneids = tr[[nch]]$clone_df[tr[[nch]]$clone_df$subject==sam,]
    databasecloneids = tr[[nch]]$clone_df[tr[[nch]]$clone_df$subject!=sam,]
    
    
    # concatenate the query clones with the closest clones in the database. 
    
    for(x in samcloneids$clone_id){
      
      queryCloneid = samcloneids[samcloneids$clone_id==x,]$did
      
      if(dbtype=="paired"){
        pw_cdr3_a_aa = tr[[nch]]$pw_cdr3_a_aa[,databasecloneids$clone_id]
        
        pw_cdr3_b_aa = tr[[nch]]$pw_cdr3_b_aa[,databasecloneids$clone_id]
        
        scoresCDr3b <- pw_cdr3_b_aa[x,]
        scoresCDr3a <- pw_cdr3_a_aa[x,]
        combinedScore <- scoresCDr3b + scoresCDr3a
        
        closeScoresidx <- which((scoresCDr3a <= 12 | scoresCDr3b <= 12) & combinedScore <= 24) # may be should use or
        
        dbHitClones = databasecloneids[databasecloneids$clone_id %in% databasecloneids$clone_id[closeScoresidx],]
        
        closeResults <- sampTCRdata_andDatabase[sampTCRdata_andDatabase$did %in% as.numeric(dbHitClones$did),]
        closeResults <- closeResults[closeResults$subject != sam,]
        
        closeResults$cdr3a_tcrdist <- scoresCDr3a[closeScoresidx]
        closeResults$cdr3b_tcrdist <- scoresCDr3b[closeScoresidx]
        closeResults$combinedPaired_tcrdist <- combinedScore[closeScoresidx]
        
        
        closeResultsDetailed <- VDJdblatest_db[match(closeResults$did,VDJdblatest_db$did),]
        closeResultsDetailed$cdr3a_tcrdist <- closeResults$cdr3a_tcrdist
        closeResultsDetailed$cdr3b_tcrdist <- closeResults$cdr3b_tcrdist
        closeResultsDetailed$combinedPaired_tcrdist <- closeResults$combinedPaired_tcrdist
        
      }else{
        #pw_cdr3_a_aa = tr[[nch]]$pw_cdr3_a_aa[,databasecloneids$clone_id]
        
        pw_cdr3_b_aa = tr[[nch]]$pw_cdr3_b_aa[,databasecloneids$clone_id]
        
        scoresCDr3b <- pw_cdr3_b_aa[x,]
        #scoresCDr3a <- pw_cdr3_a_aa[x,]
        #combinedScore <- scoresCDr3b + scoresCDr3a
        
        closeScoresidx <- which(scoresCDr3b <= 12) 
        
        dbHitClones = databasecloneids[databasecloneids$clone_id %in% databasecloneids$clone_id[closeScoresidx],]
        
        closeResults <- sampTCRdata_andDatabase[sampTCRdata_andDatabase$did %in% as.numeric(dbHitClones$did),]
        closeResults <- closeResults[closeResults$subject != sam,]
        closeResults$cdr3b_tcrdist <- scoresCDr3b[closeScoresidx]
        closeResultsDetailed <- VDJdblatest_db[match(closeResults$did,VDJdblatest_db$did),]
        closeResultsDetailed$cdr3b_tcrdist <- closeResults$cdr3b_tcrdist

      }
      
     
      #closeResultsDetailedEBV <- closeResultsDetailed[closeResultsDetailed$Pathology=="Epstein Barr virus (EBV)",]
      nEBV <- sum(closeResultsDetailed$Pathology=="Epstein Barr virus (EBV)")
      nEBVForClone <- c(nEBVForClone,nEBV)
      
      queryDetails = sampTCRdata[sampTCRdata$did==queryCloneid,]
      
      if(nrow(closeResultsDetailed) > 0){
        closeResultsDetailedComplete = closeResultsDetailed
        closeResultsDetailedComplete$Studysubject = queryDetails$subject  
        closeResultsDetailedComplete$barcode =sampTCRs[queryDetails$did,]$barcode
        closeResultsDetailedComplete$TCR = paste(queryDetails$v_a_gene,queryDetails$j_a_gene,queryDetails$cdr3_a_aa,
                                                 queryDetails$v_b_gene,queryDetails$j_b_gene,queryDetails$cdr3_b_aa,sep="_")
        
        res2 = rbind(res2,closeResultsDetailedComplete)
      }
      
      
      
    }
    
    
  }
  
  numberOfEBVhits = c(numberOfEBVhits,sum(nEBVForClone))
  
  
  write.csv(res2,file=paste0(toWrite,sam,"_",dbtype,"_VDJdb_scan.csv"))
  
}


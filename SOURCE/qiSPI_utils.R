### qiSPI ###
# description:  utils for mapping and PSM filtering
# input:        -
# output:       -
# author:       JL, HPR, GS

library(plyr)
library(dplyr)
library(stringr)

# ----- get positions -----

getPositions <- function(seq,substrate){
  
  
  
  #########################
  # PCP
  #########################
  
  
  l = nchar(seq)
  
  k = which((grepl(pattern=seq,x=substrate)==TRUE))
  if(length(k)>0){
    
    pcp = numeric()
    
    for(j in 1:length(k)){
      a = substrate
      x = strsplit(a,split=seq)[[1]]
      nn = nchar(x)
      n1 = rep(NA,(length(nn)-1))
      n2 = rep(NA,(length(nn)-1))
      for(r in 1:(length(x)-1)){
        n1[r] = sum(nn[1:r])+(r-1)*nchar(seq)+1
        n2[r] = n1[r]+nchar(seq)-1
      }
      pcp = rbind(pcp,cbind(n1,n2))
    }
    return(pcp)
  }
  
  
  
  #########################
  # PSP
  #########################
  
  
  
  
  ll = nchar(seq)
  
  
  pept = unlist(seq)
  N = nchar(seq)
  
  # split peptides to P matrix
  P = strsplit(pept,split="")[[1]]
  
  # get permutations of length N
  x = c(1:N)
  y = c(1:N)
  z = as.vector(outer(x,y,paste,sep="_"))
  q = matrix(NA,length(z),2)
  for(i in 1:length(z)){
    q[i,] = as.numeric(strsplit(z[i],split="_")[[1]])
  }
  
  qs = apply(q,1,sum)
  k = which(qs==N)
  q = q[k,]
  
  # loop over all peptides
  res2 <- list()
  res1 <- list()
  
  psp <- list()
  
  psp <- list()
  res1 <- list()
  res2 <- list()
  
  # generate all strings for searches
  S = matrix(NA,dim(q)[1],2)
  for(i in 1:dim(q)[1]){
    S[i,1] = paste(P[1:q[i,1]],sep="",collapse="")
    S[i,2] = paste(P[(q[i,1]+1):N],sep="",collapse="")
  }
  
  # search each entry in prot for the two corresponding fragments and extract acc and positions
  
  for(i in 1:dim(S)[1]){
    
    psp[[i]] <- list()
    res1[[i]] = which((grepl(pattern=S[i,1],x=substrate)==TRUE))
    res2[[i]] = which((grepl(pattern=S[i,2],x=substrate)==TRUE))
    
    kk = which(res1[[i]]%in%res2[[i]])
    k = res1[[i]][kk]
    if(length(k)>0){
      
      for(j in 1:length(k)){
        
        a = substrate
        
        
        x = strsplit(a,split=S[i,1])[[1]]
        nn = nchar(x)
        n1 = rep(NA,(length(nn)-1))
        n2 = rep(NA,(length(nn)-1))
        for(r in 1:(length(x)-1)){
          n1[r] = sum(nn[1:r])+(r-1)*nchar(S[i,1])+1
          n2[r] = n1[r]+nchar(S[i,1])-1
        }
        #check if substrate Cterm==S[i,1]
        len = nchar(S[i,1])
        y = paste(strsplit(a,split="")[[1]][(nchar(a)-len+1):nchar(a)],collapse="")
        if(S[i,1]==y){
          n1 = c(n1,nchar(a)-len+1)
          n2 = c(n2,nchar(a))
        }
        tmp = unique(apply(cbind(n1,n2),1,paste,collapse="_"))
        tmp2 = matrix(as.numeric(unlist(strsplit(tmp,split="_"))),length(tmp),2,byrow=TRUE)
        n1 = tmp2[,1]
        n2 = tmp2[,2]
        
        x = strsplit(a,split=S[i,2])[[1]]
        nn = nchar(x)
        n3 = rep(NA,(length(nn)-1))
        n4 = rep(NA,(length(nn)-1))
        for(r in 1:(length(x)-1)){
          n3[r] = sum(nn[1:r])+(r-1)*nchar(S[i,2])+1
          n4[r] = n3[r]+nchar(S[i,2])-1
        }
        #check if substrate Cterm==S[i,2]
        len = nchar(S[i,2])
        y = paste(strsplit(a,split="")[[1]][(nchar(a)-len+1):nchar(a)],collapse="")
        if(S[i,2]==y){
          n3 = c(n3,nchar(a)-len+1)
          n4 = c(n4,nchar(a))
        }
        tmp = unique(apply(cbind(n3,n4),1,paste,collapse="_"))
        tmp2 = matrix(as.numeric(unlist(strsplit(tmp,split="_"))),length(tmp),2,byrow=TRUE)
        n3 = tmp2[,1]
        n4 = tmp2[,2]
        
        # get all internal combinations and keep only those with intv<=25
        
        z = as.vector(outer(n2,n3,paste,sep="_"))
        y = matrix(NA,length(z),2)
        for(zz in 1:length(z)){
          y[zz,] = as.numeric(strsplit(z[zz],split="_")[[1]])
        }
        intv = y[,2]-y[,1]-1
        x = which(intv<0)
        if(length(x)>0){ intv[x] = y[x,1]-y[x,2]+1-nchar(S[i,1])-nchar(S[i,2]) }
        x = which(intv<0)
        #    if(length(x)>0){ intv[x] = 1000 }
        
        select = which(intv<=5000)
        
        nnn = length(select)
        if(nnn>0){
          psp[[i]][[j]] = matrix(NA,nnn,5)
          
          for(j2 in 1:nnn){
            
            psp[[i]][[j]][j2,] = c(pept,y[select[j2],1]-nchar(S[i,1])+1,y[select[j2],1],y[select[j2],2],y[select[j2],2]+nchar(S[i,2])-1)
          }
        }
        
      }
      
    }
    
    
  }
  
  # unlist results and return as unique matrix with all possible explanations as rows
  x = unlist(psp)
  #  psp = matrix(x,length(x)/5,5,byrow=FALSE)
  
  res = numeric()
  for(i in 1:length(psp)){
    if(length(psp[[i]])>0){
      for(j in 1:length(psp[[i]])){
        res = rbind(res,psp[[i]][[j]])
      }
    }
  }
  
  
  # print(res)
  
  return(res)
  
  
}

# ----- re-map peptides -----
# re-map peptides
# substrateSeq, productType, spliceType, pepSeq

mapping = function(DB) {
  
  d = DB %>%
    select(substrateSeq, productType, spliceType, pepSeq) %>%
    distinct()
  
  pb = txtProgressBar(min = 0, max = dim(d)[1], style = 3)
  
  for(i in 1:dim(d)[1]){
    setTxtProgressBar(pb, i)
    
    if(!(d$productType[i]=="CONT" | d$productType[i]=="CONT_synError")){
      s = gsub("I","L",as.vector(d$pepSeq[i]))
      substrate = gsub("I","L",as.vector(d$substrateSeq[i]))
      x = getPositions(s,substrate)
      
      
      #PCP
      if(dim(x)[2]==2){
        d$positions[i] = paste(apply(x,1,paste,collapse="_"),collapse=";")
      }
      
      
      #PSP
      if(dim(x)[2]>2){
        # print(i)
        if(dim(x)[2]==5 & dim(x)[1]>1){
          d$positions[i] = paste(apply(x[,-1],1,paste,collapse="_"),collapse=";")
        }
        if(dim(x)[2]==5 & dim(x)[1]==1){
          d$positions[i] = paste(x[,-1],collapse="_")
        }
        
        types = rep("cis",dim(x)[1])
        
        intv = as.numeric(x[,4])-as.numeric(x[,3])
        k = which(intv<=0)
        if(length(k)>0){
          types[k] = "revCis"
          
          k2 = which(as.numeric(x[k,2])<=as.numeric(x[k,5]) & as.numeric(x[k,3])>=as.numeric(x[k,4]))
          if(length(k2)>0){
            types[k[k2]] = "trans"
          }
          
        }
        
        d$spliceType[i] = paste(types,collapse=";")
      }
    }
    
  }
  
  DB$productType = NULL
  DB$spliceType = NULL
  DB$positions = NULL
  
  
  DB = left_join(DB, d) %>%
    as.data.frame()
  
  # re-assign the product type
  pos = strsplit(DB$positions, "_")
  wr = which(sapply(pos, length) == 2 & str_detect(DB$productType, "PSP"))
  if (length(wr) > 0) {
    DB$productType[wr] = str_replace(DB$productType[wr], "PSP", "PCP")
  }
  
  xr = which(sapply(pos, length) == 4 & str_detect(DB$productType, "PCP"))
  if (length(xr) > 0) {
    DB$productType[xr] = str_replace(DB$productType[xr], "PCP", "PSP")
  }
  
  return(DB)
}


# ----- PSM filtering -----
# input all PSMs that were detected, with the same column names as in invitroSPI
# scanNum, runID, rank, pepSeq, productType, rank, ionScore, qValue, ...

PSMfiltering = function(MSDB, delta_score, ion_score, q_value) {
  
  columns = names(MSDB)
  runIDs = unique(MSDB$runID)
  
  allPSMs = list()
  allPSMs_Delta = list()
  
  pb = txtProgressBar(min = 0, max = length(runIDs), style = 3)
  for (i in 1:length(runIDs)) {
    
    setTxtProgressBar(pb, i)
    
    selected <- data.frame(matrix(ncol = length(columns), nrow = 0))
    colnames(selected) = columns
    
    selectedDeltaRecord = data.frame(matrix(ncol = length(columns), nrow = 0))
    colnames(selectedDeltaRecord) <- columns
    
    runIDTable = MSDB[MSDB$runID == runIDs[i],]
    scanNum = as.character(runIDTable$scanNum) %>% unique()
    
    
    # Iterate through scanNum and sort by rank
    for (k in 1:length(scanNum)) {
      
      PSPcandidatesIndex = NULL
      PCPcandidatesIndex = NULL
      
      deltaRecord = data.frame(matrix(ncol = length(columns), nrow = 0))
      colnames(deltaRecord ) = columns
      
      scanNumTable = runIDTable[runIDTable$scanNum == scanNum[k],]
      scanNumTable = scanNumTable[order(scanNumTable$rank),]
      
      # # How many entries in filteredScans are rank 1?
      scanNumTable = unique(scanNumTable)
      filteredScans = scanNumTable
      filteredScans = filteredScans[order(filteredScans$rank),]
      
      # print(filteredScans)
      
      # extract the top score
      topScore <- filteredScans[1, "ionScore"]
      
      # extract scans with rank 1
      allRanks = table(filteredScans$rank)
      NoTopRankedScans <- allRanks[names(allRanks) == 1]
      
      # it can somehow happen that there is no top-ranked peptide...
      if (length(NoTopRankedScans) == 0) {
        filteredScans <- filteredScans[0,]
        paste0("!!! NO RANK 1 PEPTIDES FOUND IN SCAN ", scanNum[k], ", RUN-ID ", runIDs[i], " !!!") %>%
          print()
        
      } else if (NoTopRankedScans > 1) {
        
        topIndices = which(filteredScans$rank == 1)
        NoTopRankedPeps = gsub("I","L",filteredScans$pepSeq[topIndices]) %>%
          unique() %>%
          length()
        
        
        PCPindex = which(filteredScans$productType[filteredScans$rank == 1]=="PCP")
        PSPindex = which(filteredScans$productType[filteredScans$rank == 1]=="PSP")
        
        # if more than one peptide is rank 1
        # keep scan only if there is just a single PCP, otherwise discard
        if (NoTopRankedPeps > 1) {
          
          if (length(PCPindex) >= 1) {
            
            if (length(unique(gsub("I","L",filteredScans$pepSeq[PCPindex]))) == 1) {
              filteredScans = filteredScans[PCPindex[1],]
            } else {
              filteredScans = filteredScans[0,]
            }
            
          } else {
            
            if (length(PSPindex) > 1) {
              index = which(filteredScans$productType[which(filteredScans$rank == 1)]=="PSP")
              PSPcandidatescores <- filteredScans[index, "ionScore"]
              ind = which(1-PSPcandidatescores/topScore <= delta_score)
              deltaRecord = rbind(deltaRecord,filteredScans[index[ind],])
            }
            
            filteredScans = filteredScans[0,]
            
          }
          
        } else {
          filteredScans = filteredScans[1,]
        } 
        
      }
      # if there is only one scan on rank 1, there is also only one peptide
      
      # if there is only a single scan on rank 1
      # (nrow(filteredScans) will be larger than 1 since the previous if condition did not apply)
      # if the 1st rank is a PCP, assign it as such
      else if (NoTopRankedScans == 1 & filteredScans[1, "productType"] == "PCP") {
        
        filteredScans <- filteredScans[1,]
        
        # If rank 1 entry is PSP check if a likely PCP is present in the lower ranks (PCPcandidates)
        # also check if there are lower-ranked PSPs
      } else if (NoTopRankedScans == 1 & filteredScans[1, "productType"] == "PSP") {
        
        # lower-ranked PCPs and PSPs
        PCPcandidatesIndex_temp <- which(as.numeric(filteredScans$rank) > 1 & as.character(filteredScans$productType) == "PCP")
        PSPcandidatesIndex_temp <- which(as.numeric(filteredScans$rank) > 1 & as.character(filteredScans$productType) == "PSP")
        
        # if there are no PSPs or PCPs with a lower rank, assign the PSP
        if(length(PCPcandidatesIndex_temp)==0 & length(PSPcandidatesIndex_temp)==0){
          filteredScans <- filteredScans[1,]
        }
        
        # if there are PSPs with a lower rank --> get highest-ranked PSM
        if (length(PSPcandidatesIndex_temp)>0) {
          #keeps only top PSP candidate
          PSPcandidatesIndex = which(as.numeric(filteredScans$rank)==min(as.numeric(filteredScans$rank)[PSPcandidatesIndex_temp]) & as.character(filteredScans$productType) == "PSP")[1]
          # keeps all PSP candidates
          PSPcandidatesIndexAll = PSPcandidatesIndex_temp
        } else {
          PSPcandidatesIndex = NULL
        }
        
        # if there are PCPs with a lower rank --> get highest-ranked PSM
        if (length(PCPcandidatesIndex_temp)>0) {
          PCPcandidatesIndex = which(as.numeric(filteredScans$rank)==min(as.numeric(filteredScans$rank)[PCPcandidatesIndex_temp]) & as.character(filteredScans$productType) == "PCP")
        } else {
          PCPcandidatesIndex = NULL
        }
        
        # for lower-ranked PCPs --> calculate Delta score
        # Keep top ranked PSP if deltascore > 0.3 and remove rest
        if (length(PCPcandidatesIndex) > 0) {
          
          
          PCPcandidatescore <- filteredScans[PCPcandidatesIndex[1], "ionScore"]
          # if the Delta score between the top-ranked PSP and a lower-ranked PCP
          # is larger than the threshold
          if (1 - PCPcandidatescore / topScore > delta_score) {
            
            # if there are any lower-ranked PSPs
            # Calculate delta_score scores also for lower-ranked PSPs
            if (length(PSPcandidatesIndex)>0) {
              PSPcandidatescore <- filteredScans[PSPcandidatesIndex, "ionScore"]
              
              # lower-ranked PSP, but with very low ion score --> assign top-ranked PSP
              if (1 - PSPcandidatescore / topScore > delta_score) {
                filteredScans <- filteredScans [1,]
                
                # discard spectrum bc there are several PSP candidates with similar
                # ion score
              } else if (1 - PSPcandidatescore / topScore <= delta_score) {
                
                PSPcandidatescores <- filteredScans[PSPcandidatesIndexAll, "ionScore"]
                ind = which(1-PSPcandidatescores/topScore <= delta_score)
                deltaRecord = rbind(deltaRecord,filteredScans[PSPcandidatesIndexAll[ind],])
                
                filteredScans <- filteredScans[0,]
              }
              
              # if there are no lower-ranked PSPs --> assign the top-ranked PSP
            } else if (length(PSPcandidatesIndex)==0){
              filteredScans <- filteredScans[1,]
            }
            
            # if the Delta score is too small:
            # Take candidate PCP as new top ranked and remove rest
          } else if (1 - PCPcandidatescore / topScore <= delta_score) {
            if (length(PCPcandidatesIndex) == 1) {
              filteredScans <- filteredScans[PCPcandidatesIndex,]
            } else {
              
              # more than one lower-ranked PCP that survived Delta filter
              if (length(unique(gsub("I","L",filteredScans$pepSeq[PCPcandidatesIndex]))) == 1) {
                filteredScans <- filteredScans[PCPcandidatesIndex[1],]
              } else {
                filteredScans = filteredScans[0,]
              }
              
            }
            
          }
          
          # end if there are both PCPs and PSPs in lower ranks
          # if there are lower ranked PSPs but no PCPs: continue
          
        } else if (length(PSPcandidatesIndex) > 0 & length(PCPcandidatesIndex)==0){
          PSPcandidatescore <- filteredScans[PSPcandidatesIndex, "ionScore"]
          
          # keep the candidate PSP if the Delta score is large enough
          if (1 - PSPcandidatescore / topScore > delta_score) {
            filteredScans <- filteredScans [1,]
            
            # if the Delta score is too small --> discard entire scan
          } else if (1 - PSPcandidatescore / topScore <= delta_score) {
            
            PSPcandidatescores <- filteredScans[PSPcandidatesIndexAll, "ionScore"]
            ind = which(1-PSPcandidatescores/topScore <= delta_score)
            deltaRecord = rbind(deltaRecord,filteredScans[PSPcandidatesIndexAll[ind],])
            
            filteredScans <- filteredScans[0,]
          }
          
        }
      }
      
      # apply filter for q-value and ion score
      if(nrow(filteredScans) > 0){
        if(as.double(filteredScans$qValue) >= q_value | as.double(filteredScans$ionScore <= ion_score)){
          filteredScans <- filteredScans[0,]
        }else{
          selected <- rbind(filteredScans, selected)
        }
      }
      
      if(nrow(deltaRecord)> 0){
        k = which(as.double(deltaRecord$qValue) < q_value & as.double(deltaRecord$ionScore > ion_score))
        deltaRecord <- deltaRecord[k,]
        
        if(nrow(deltaRecord)> 0){
          selectedDeltaRecord <- rbind(deltaRecord,selectedDeltaRecord)
        }
      }
      
    }
    
    allPSMs[[i]] = selected
    allPSMs_Delta[[i]] = selectedDeltaRecord
    
  }
  
  allPSMs = plyr::ldply(allPSMs)
  allPSMs_Delta = plyr::ldply(allPSMs_Delta)
  
  print("")
  paste0("number of PSMs: ", nrow(allPSMs)) %>%
    print()
  
  return(list(allPSMs = allPSMs,
              allPSMs_Delta = allPSMs_Delta))
}



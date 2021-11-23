### iqSPI ###
# description:  reassign types of quantitation queries
# input:        sample list, psps and pcps, quantitation results, num queries
# output:       re-assigned quantitation results
# author:       JL, modified by HPR

library(seqinr)
library(dplyr)

print("---------------------------------")
print("4) RE-ASSIGN QUANTITATION QUERIES")
print("---------------------------------")

protein_name = snakemake@params[["protein_name"]]
PROTEINSEQUENCES_PATH = "INPUT/sequences/"


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name %in% protein_name, ]

load(snakemake@input[["numQueries"]])


### MAIN PART ###

prots = sample_list$protein_name %>%
  unique()

for (p in prots) {
  paste0("re-assigning queries of ", p) %>%
    print()
  
  cnt = sample_list[sample_list$protein_name == p,]
  cntNumQueries = numQueries[p] %>%
    unlist()
  
  for(i in 1:nrow(cnt)){
    
    indir = gsub(".raw","",cnt$raw_file[i])
    
    #----- load search results PCP and PSP -----
    load(file=paste0("OUTPUT/tmp/",p,"/",indir,"/pcp.RData"))
    load(file=paste0("OUTPUT/tmp/",p,"/",indir,"/psp.RData"))
    
    pcpQuery = as.numeric(as.vector(pcp$pep_query))
    pspQuery = as.numeric(as.vector(psp$pep_query))
    
    pcpSeq = as.vector(pcp$pep_seq)
    pspSeq = as.vector(psp$pep_seq)
    
    pcpACC = as.vector(pcp$prot_acc)
    pspACC = as.vector(psp$prot_acc)
    
    
    #----- load quantitation results -----
    quant = read.csv(paste0("OUTPUT/tmp/",p,"/",indir,"/quantitationResults.csv"),
                     stringsAsFactors = F)
    
    if(i>1){
      quantQuery = as.numeric(as.vector(quant$Query))-sum(as.numeric(cntNumQueries[1:(i-1)]))
    }
    if(i==1){
      quantQuery = as.numeric(as.vector(quant$Query))
    }
    
    #----- reassign -----
    
    quantAccession = rep(NA,dim(quant)[1])
    quantSequence = rep(NA,dim(quant)[1])
    
    for(j in 1:dim(quant)[1]){
      
      k1 = which(pcpQuery==quantQuery[j])
      k2 = which(pspQuery==quantQuery[j])
      
      if(length(k1)==1){
        quantAccession[j] = pcpACC[k1]
        quantSequence[j] = pcpSeq[k1]
      }
      
      if(length(k2)==1){
        quantAccession[j] = pspACC[k2]
        quantSequence[j] = pspSeq[k2]
      }
      
      
    }
    
    sequence = quantSequence
    quant$Accession = quantAccession
    quant = cbind(quant,sequence)
    
    write.csv(quant,
              file=paste0("OUTPUT/tmp/",p,"/",indir,"/reassigned.csv"),
              row.names = F)
    
  }
  
}


### OUTPUT ###
write("queries in quantitation results were reassigned",
      file = "OUTPUT/tmp/reassign_queries.txt")


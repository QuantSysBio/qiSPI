### iqSPI ###
# description:  reassign types of quantitation queries
# input:        sample list, re-assigned quantitation results, protein sequences
# output:       mapping (merged over replicates)
# author:       JL, modified by HPR

library(seqinr)
library(dplyr)

print("-------------------------------")
print("5) MERGE RE-ASSIGNED REPLICATES")
print("-------------------------------")

protein_name = snakemake@params[["protein_name"]]
PROTEINSEQUENCES_PATH = "INPUT/sequences/"

PCPthresh = snakemake@params[["PCPthresh"]]
PSPthresh = snakemake@params[["PSPthresh"]]


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name %in% protein_name, ]


### MAIN PART ###
prots = sample_list$protein_name %>%
    unique()

for (p in prots) {
    
    paste0("merging all replicates of ", p) %>%
        print()
    
    cnt = sample_list[sample_list$protein_name == p,]
    FF = gsub(".raw", "", cnt$raw_file)
    
    #----- load reassigned tables per replicate -----
    rep = list()
    for(fileInd in 1:length(FF)){
        
        tmp = read.csv(paste0("OUTPUT/tmp/",p,"/",FF[fileInd],"/reassigned.csv"),
                       stringsAsFactors = F)
        k = which(is.na(tmp$Accession)==FALSE)
        rep[[fileInd]] = tmp[k,]
        
    }
    
    numReplicates = length(FF)
    
    #----- combine replicate tables: protHit, acc1, seq1 -----
    protHit = list()
    acc = list()
    seq = list()
    charge = list()
    
    for(i in 1:numReplicates){
        protHit[[i]] = as.numeric(as.vector(rep[[i]]$protHit))
        acc[[i]] = as.vector(rep[[i]]$Accession)
        seq[[i]] = as.vector(rep[[i]]$sequence)
        charge[[i]] = as.numeric(as.vector(rep[[i]]$z))
    }
    
    allprotHit = numeric()
    for(i in 1:numReplicates){
        allprotHit = c(allprotHit,paste(protHit[[i]],seq[[i]],sep="_"))
    }
    
    protHitunique = sort(unique(allprotHit))
    
    x = unlist(lapply(protHitunique,strsplit,split="_"))
    protHitunique = as.numeric(x[seq(1,length(x),2)])
    seqHitunique = x[seq(2,length(x),2)]
    
    #----- get sequences for unique prot_hits -----
    if (grepl(".fasta", cnt$substrateSeq[1])) {
        protInfile = paste0(PROTEINSEQUENCES_PATH, "/", cnt$substrateSeq[1])
        proteome = read.fasta(file=protInfile,seqtype = "AA", as.string=TRUE)[1] %>% unlist()
    } else {
        proteome = cnt$substrateSeq[1]
    }
    
    #----- combine all tables -----
    mapping = matrix(NA,length(protHitunique),3)
    
    for(i in 1:length(protHitunique)){
        
        mapping[i,1] = protHitunique[i]
        mapping[i,3] = seqHitunique[i]
        
        x = list()
        y = list()
        for(j in 1:numReplicates){
            k = which(protHit[[j]]==protHitunique[i])
            if(length(k)>0){
                x[[j]] = acc[[j]][k]
                y[[j]] = seq[[j]][k]
            }
        }
        
        t1 = table(unlist(x))
        t2 = table(unlist(y))
        
        if(length(t2)>0){
            r = t2/sum(t2)*100
            n = rownames(t2)
            
            # PCP in n?
            ind = numeric()
            for(j in 1:length(n)){
                k = grep(n[j],proteome)
                if(length(k)>0){ ind = c(ind,j) }
            }
            
            if(length(ind)>0){
                ind = which(r[ind]==max(r[ind]))[1]
                if(r[ind]>PCPthresh){
                    mapping[i,2] = n[ind]
                }
            }
            
            # no PCP in n
            if(length(ind)==0){
                ind = which(r==max(r))[1]
                if(r[ind]>PSPthresh){
                    mapping[i,2] = n[ind]
                }
                
            }
        }
        
    }
    
    # ----- save mappings -----
    outdir = paste0("OUTPUT/",p,"/")
    if(!dir.exists(outdir)) {
        dir.create(outdir)
    }
    
    save(mapping,
         file=paste0(outdir,"mapping.RData"))
    
}


### OUTPUT ###
write("replicates were merged",
      file = "OUTPUT/tmp/merge_replicates.txt")


### qiSPI ###
# description:  post-processing of kinetics
# input:        sample list, quantification results per biological replicate,
#               protein sequence, psp.RData, pcp.RData
# output:       processed kinetics
# author:       JL, modified by HPR


library(dplyr)
library(seqinr)
library(stringr)

print("------------------------------")
print("8) POST-PROCESSING OF KINETICS")
print("------------------------------")

protein_name = snakemake@params[["protein_name"]]
PROTEINSEQUENCES_PATH = "INPUT/sequences/"


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name %in% protein_name, ]


### MAIN PART ###

prots = sample_list$protein_name %>%
    unique()

for (p in prots) {
    
    paste0("processing kinetics of ", p) %>%
        print()
    
    cnt = sample_list[sample_list$protein_name == p,]
    t = cnt$digestTime %>%
        unique() %>%
        as.numeric() %>%
        sort(decreasing = F)
    
    
    # ----- load splitted intensities and protein sequence -----
    load(paste0("OUTPUT/",p,"/quantity_rep.RData"))
    rep = lapply(rep, function(x) { x[,-1] })
    
    si = cnt$substrateSeq[1]
    if (grepl(".fasta", si)) {
        protInfile = paste0(PROTEINSEQUENCES_PATH, "/", si)
        proteome = read.fasta(file=protInfile,seqtype = "AA", as.string=TRUE)[1] %>% unlist()
    } else {
        proteome = si
    }
    
    # ----- split into technical replicates -----
    counter = 1
    reps = list()
    c = which(names(rep[[1]]) == "z")
    
    for (r in 1:length(rep)) {
        R = rep[[r]]
        tr = names(R)[grepl("tp_", names(R))]
        
        techReps = str_split(tr, pattern = coll("."), simplify = T)[,1] %>%
            table()
        techReps = techReps[1]
        
        for (tt in 1:techReps) {
            cidx = ((c+1):ncol(R))-c
            cx = seq(tt,length(cidx),techReps)
            reps[[counter]] = R[,c(1:c, c+cx)]
            
            counter = counter + 1
        }
    }
    
    names(reps) = c(1:(counter-1))
    numReps = length(reps)
    paste0("number of replicates: ", numReps) %>%
        print()
    
    # ----- aggregate all charges and PTMs by summation -----
    su = unique(as.vector(rep[[1]]$sequence))
    
    results = list()
    for(i in 1:numReps){
        results[[i]] = matrix(NA,length(su),(length(t)+1))
    }
    
    for(i in 1:length(su)){
        for(j in 1:numReps){
            
            k = which(as.vector(reps[[j]]$sequence)==su[i])
            results[[j]][,1] = su
            if(length(k)>1){
                results[[j]][i,2:(length(t)+1)] = as.numeric(as.vector(apply(reps[[j]][k,(c+1):(length(t)+c)],2,function(x) {
                    sum(as.numeric(x))
                })))
            }
            if(length(k)==1){
                results[[j]][i,2:(length(t)+1)] = as.numeric(as.vector(reps[[j]][k,(c+1):(length(t)+c)]))
            }
        }
    }
    
    
    paste0("size of results: ", length(results)) %>%
        print()
    
    # ----- get all peptides detected at t=0 -----
    ctrlInd = which(cnt$digestTime == 0 | cnt$digestTime == "CTRL")
    
    # only remove if control measurements are present
    if (length(ctrlInd) > 0) {
        
        FF = gsub(".raw","",cnt$raw_file)
        
        allPCP = c()
        allPSP = c()
        for(i in 1:length(ctrlInd)){
            load(file=paste0("OUTPUT/tmp/",p,"/",FF[ctrlInd[i]],"/psp.RData"))
            allPSP = c(allPSP,as.vector(psp$pepSeq))
        }
        ctrlPSP = unique(gsub("I","L",unique(allPSP)))
        
        # ----- remove all control peptides as such -----
        k2 = which(gsub("I","L",results[[1]][,1]) %in% ctrlPSP)
        k = k2
        if(length(k)>0){
            for(i in 1:numReps){
                results[[i]] = results[[i]][-k,]
            }
        }
        
    }
    
    # ----- assign types: PCP vs PSP -----
    pcp = c()
    psp = c()
    for(i in 1:dim(results[[1]])[1]){
        k = grep(gsub("I","L",results[[1]][i,1]),gsub("I","L",proteome))
        if(length(k)>0){
            pcp = c(pcp,i)
        }
        if(length(k)==0){
            psp = c(psp,i)
        }
    }
    
    types = rep(NA,dim(results[[1]])[1])
    types[pcp] = "pcp"
    types[psp] = "psp"
    
    
    # ----- remove all PSPs that have a precursor as ctrlPSP -----
    if (length(ctrlInd) > 0) {
        remove = c()
        for(i in psp){
            s = gsub("I","L",results[[1]][i,1])
            k = grep(s,gsub("I","L",ctrlPSP))
            if(length(k)>0){
                remove = c(remove,i)
            }
        }
        
        if(length(remove)>0){
            
            types = types[-remove]
            for(i in 1:numReps){
                results[[i]] = results[[i]][-remove,]
            }
            
        }
        
    }
    
    # ------ format and save results -----
    cnames = names(reps[[1]])[str_detect(names(reps[[1]]), "tp_")]
    cnames = str_split(cnames, "_", simplify = T)[,2] %>%
        str_replace_all(pattern = ",", replacement = ".") %>%
        as.numeric()
    cord = order(cnames, decreasing = F)
    
    for(i in 1:numReps){
        results[[i]] = data.frame(cbind(results[[i]],types))
        results[[i]] = results[[i]][, c(1,cord+1,ncol(results[[i]]))]
        
        colnames(results[[i]]) = c("sequence",t,"types")
    }
    
    save(results,file=paste0("OUTPUT/",p,"/results.RData"))
    
}


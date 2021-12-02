## iqSPI ###
# description:  assign product types to filtered means/results, plot final kinetics
# input:        filtered results/means, protein sequences
# output:       final filtered results/means
# author:       JL, modified by HPR

library(dplyr)
library(seqinr)
source("SOURCE/loadFunctions.r")

print("----------------------------------------------------------")
print("10) COMPLETE TYPE ASSIGNMENT, FINAL KINETICS AND MERGED DB")
print("----------------------------------------------------------")

protein_name = snakemake@params[["protein_name"]]
rm_tmp = snakemake@params[["rm_tmp"]]
PROTEINSEQUENCES_PATH = "INPUT/sequences/"


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name %in% protein_name, ]


### MAIN PART ###
prots = sample_list$protein_name %>%
    unique()

allKinetics = list()
counter = 1
for (p in prots) {
    
    paste0("processing kinetics of ", p) %>%
        print()
    
    cnt = sample_list[sample_list$protein_name == p,]
    t = cnt$digestTime %>%
        unique() %>%
        as.numeric() %>%
        sort(decreasing = F)
    
    # ----- load filtered results/means -----
    load(file=paste0("OUTPUT/",p,"/filteredResults.RData"))
    load(file=paste0("OUTPUT/",p,"/filteredMeans.RData"))
    
    M = filteredMeans
    
    si = cnt$substrateSeq[1]
    if (grepl(".fasta", si)) {
        protInfile = paste0(PROTEINSEQUENCES_PATH, "/", si)
        proteome = read.fasta(file=protInfile,seqtype = "AA", as.string=TRUE)[1] %>% unlist()
    } else {
        proteome = si
    }
    proteome = gsub("I", "L", proteome)
    
    # ----- replace all Is with Ls to find alternative sequences -----
    seqSub <- gsub("I", "L", as.vector(filteredResults[[1]]$sequence))
    
    # ----- get assignment ----
    positions = rep(NA,length(seqSub))
    spliceType = rep("non-spliced",length(seqSub))
    
    for (i in 1:length(seqSub)){
        
        
        s = gsub("I","L",as.vector(seqSub[i]))
        substrate = gsub("I","L",as.vector(proteome))
        x = getPositions(s,substrate)
        
        #PCP
        if(dim(x)[2]==2){
            positions[i] = paste(apply(x,1,paste,collapse="_"),collapse=";")
        }
        
        #PSP
        if(dim(x)[2]>2){
            if(dim(x)[2]==5 & dim(x)[1]>1){
                positions[i] = paste(apply(x[,-1],1,paste,collapse="_"),collapse=";")
            }
            if(dim(x)[2]==5 & dim(x)[1]==1){
                positions[i] = paste(x[,-1],collapse="_")
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
            spliceType[i] = paste(types,collapse=";")
            
        }
    }
    
    # ----- format results -----
    for(i in 1:length(filteredResults)){
        filteredResults[[i]] = cbind(filteredResults[[i]],positions,spliceType)
    }
    
    # ----- plot final kinetics -----
    pdf(paste0("OUTPUT/",p,"/PLOTS/finalKinetics.pdf"), width=10, height=10)
    par(mfrow = c(4,5))
    
    myColors = c("red","pink","blue","lightblue")
    for(i in 1:dim(M[[1]])[1]){
        
        maxi = 0
        mini = 10**20
        for(j in 1:length(M)){
            maxi = max(c(maxi,as.numeric(as.vector(M[[j]][i,]))))
            mini = min(c(mini,as.numeric(as.vector(M[[j]][i,]))))
        }
        
        plot(1,1,
             type="b",
             col="white",
             xlim=c(0,max(t)),axes=FALSE,
             main=filteredResults[[1]]$positions[i],
             ylim=c(mini,maxi),xlab="time",ylab="Intensity")
        axis(1)
        axis(2)
        
        for(j in 1:length(M)){
            points(t,M[[j]][i,],type="b",col=myColors[j])
        }
        
    }
    
    dev.off()
    
    # ----- format and save results -----
    for(i in 1:length(M)){
        filteredMeans[[i]] = cbind(filteredResults[[i]]$sequence,M[[i]],filteredResults[[i]]$types, positions,spliceType)
    }
    
    filteredMeans = lapply(filteredMeans, function(x) {
        x = as.data.frame(x)
        names(x) = c("pepSeq", paste0("tp_",t), "productType", "positions", "spliceType")
        return(x)
    })
    names(filteredMeans) = c(1:length(filteredMeans))
    
    # kinetics as data frame
    Kinetics.DF = plyr::ldply(filteredMeans)
    names(Kinetics.DF)[1] = "biological_replicate"
    Kinetics.DF$substrateID = cnt$substrateID[1]
    Kinetics.DF$substrateSeq = cnt$substrateSeq[1]
    
    o = order(Kinetics.DF$substrateID,Kinetics.DF$pepSeq, Kinetics.DF$biological_replicate)
    Kinetics.DF = Kinetics.DF[o,]
    
    Kinetics.DF.pasted = Kinetics.DF
    Kinetics.DF.pasted$digestTimes = paste(t, collapse = ";")
    xx = which(grepl("tp_", names(Kinetics.DF.pasted)))
    Kinetics.DF.pasted$intensities = do.call(paste, c(Kinetics.DF.pasted[xx], sep=";"))
    
    allKinetics[[counter]] = Kinetics.DF.pasted[,-which(grepl("tp_", names(Kinetics.DF.pasted)))]
    
    # peptides only
    pepList = allKinetics[[counter]][,c("substrateID","substrateSeq","pepSeq",
                                        "productType", "spliceType", "positions",
                                        "digestTimes")] %>% unique()
    
    counter = counter + 1
    
    ### OUTPUT ###
    save(filteredResults,
         file=paste0("OUTPUT/",p,"/filteredResults_final.RData"))
    save(filteredMeans,
         file=paste0("OUTPUT/",p,"/filteredMeans_final.RData"))
    write.csv(Kinetics.DF,
              file = paste0("OUTPUT/",p,"/finalKinetics.csv"),
              row.names = F)
    write.csv(pepList,
              file = paste0("OUTPUT/",p,"/pepList.csv"),
              row.names = F)
    
}


# ---- construct kinetics DB -----
names(allKinetics) = prots
allKinetics = plyr::ldply(allKinetics)
names(allKinetics)[1] = "protein_name"

write.csv(allKinetics,
          file = unlist(snakemake@output[["KineticsDB"]]),
          row.names = F)


# ---- remove temporary directories -----

if (rm_tmp == "yes") {
    system("rm -rf OUTPUT/tmp/")
}


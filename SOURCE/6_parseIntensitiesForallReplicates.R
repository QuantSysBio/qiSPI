### qiSPI ###
# description:  parse intensities for all replicates
# input:        sample list, peptide quantification results, mapping
# output:       quantities (merged over replicates)
# author:       JL, modified by HPR

library(seqinr)
library(dplyr)

print("---------------------------------------")
print("6) PARSE INTENSITIES FOR ALL REPLICATES")
print("---------------------------------------")

protein_name = snakemake@params[["protein_name"]]
QUANTITATIONRESULTS_PATH = "INPUT/quantitation_results/"


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name %in% protein_name, ]


### MAIN PART ###
qpeps = sample_list$peptide_quantification %>%
    unique()

for (q in qpeps) {
    
    paste0("parsing all intensities of ", q) %>%
        print()
    
    
    #----- load peptide quantification results -----
    infile = paste0(QUANTITATIONRESULTS_PATH,q)
    if (grepl(".csv", infile)) {
        quant = read.csv(infile, stringsAsFactors = F)
    } else if (grepl(".txt", infile)) {
        quant = read.csv(infile, stringsAsFactors = F, sep = "\t")
    }
    
    FF = gsub(".raw","",sample_list$raw_file[sample_list$peptide_quantification == q])
    numReplicates = length(FF)
    
    
    #----- load mapping -----
    p = sample_list$protein_name[sample_list$peptide_quantification == q][1]
    load(paste0("OUTPUT/",p,"/mapping.RData"))
    
    x = colnames(quant)
    k = grep("Intensity",x)
    quant = quant[,c(1:2,7,k)]
    quant[,2] = NA
    
    
    #----- remove not assigned protHits -----
    
    k = which(duplicated(mapping[,2]) | is.na(mapping[,2]))
    if(length(k)>0){mapping = mapping[-k,]}
    
    
    #----- add correct assignment -----
    for(i in 1:dim(quant)[1]){
        
        k = which(mapping[,1]==quant$Hit[i])
        if(length(k)>0){
            quant[i,1:2] = mapping[k,1:2]
        }
        
    }
    
    k = which(is.na(quant[,2])==TRUE)
    if(length(k)>0){
        quant = quant[-k,]
    }
    
    
    #----- make quantity matrix -----
    quantity = matrix(NA,dim(quant)[1],dim(quant)[2])
    rem = numeric()
    
    for(i in 1:dim(quant)[1]){
        
        k = which((quant[,1]==quant$Hit[i])&(quant[,3]==quant$z[i]))
        
        avg = apply(quant[k,4:dim(quant)[2]],2,mean)
        quantity[i,c(4:dim(quant)[2])] = avg
        quantity[i,1] = as.vector(quant[i,1])
        quantity[i,2] = as.vector(quant[i,2])
        quantity[i,3] = as.vector(quant[i,3])
        
        if(length(k)>1){
            rem = c(rem,k[-1])
        }
        
    }
    
    rem = unique(rem)
    if(length(rem)>0){
        quantity = quantity[-rem,]
    }
    
    
    #----- rename as data frame ------
    quantity = data.frame(quantity)
    colnames(quantity)<-c("protHit","sequence","z",paste("rep_",1:numReplicates,sep=""))
    
    
    ### OUTPUT ###
    save(quantity,file=paste0("OUTPUT/",p,"/quantity.RData"))
    write.csv(quantity,
              file=paste0("OUTPUT/",p,"/quantity.csv"),
              row.names = F)
    
}


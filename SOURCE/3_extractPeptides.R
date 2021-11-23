### iqSPI ###
# description:  extract peptides from search result file
# input:        sample list (protein sequences), parsed search results
# output:       PCPs and PSPs per replicate
# author:       JL, modified by HPR

library(seqinr)
library(stringr)
library(dplyr)

print("--------------------------------------------")
print("3) EXTRACT PEPTIDES FROM SEARCH RESULT FILES")
print("--------------------------------------------")

protein_name = snakemake@params[["protein_name"]]
PROTEINSEQUENCES_PATH = "INPUT/sequences/"

KK = snakemake@params[["KK"]]
IONscore = snakemake@params[["IONscore"]]
thrDif = snakemake@params[["thrDif"]]
thrDifpcp = snakemake@params[["thrDifpcp"]]

PSPstring = "PSP_"
PCPstring = "PCP_"


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name %in% protein_name, ]

load(snakemake@input[["rawNames"]])


### MAIN PART ###
seqs = sample_list$substrateSeq %>%
    unique()

for (si in seqs) {
    
    print("--------------------------")
    paste0("LOADING: ", si) %>% print()
    print("--------------------------")
    
    if (grepl(".fasta", si)) {
        protInfile = paste0(PROTEINSEQUENCES_PATH, "/", si)
        proteome = read.fasta(file=protInfile,seqtype = "AA", as.string=TRUE)[1] %>% unlist()
    } else {
        proteome = si
    }
    
    proteome = gsub("I","L",proteome)
    
    cnt = sample_list[sample_list$substrateSeq == si,]
    
    for (c in 1:nrow(cnt)) {
        
        # ----- load search results -----
        cntRawNames = rawNames[names(rawNames) == cnt$protein_name[1]] %>%
            unlist()
        names(cntRawNames) = str_split(names(cntRawNames), coll("."), simplify = T)[,2] %>%
            as.character()
        indir = names(cntRawNames)[cntRawNames == cnt$order[c] & names(cntRawNames) %in% gsub(".raw","",cnt$raw_file)]
        print(indir)
        
        dat = read.csv(paste0("OUTPUT/tmp/",cnt$protein_name[c],"/",indir,"/searchResults.csv"),
                       stringsAsFactors = F)
        
        i = order(dat$pep_query)
        dat = dat[i,]
        
        
        # ----- remove empty queries -----
        seq = as.vector(dat$pep_seq)
        ind = which(is.na(seq)==FALSE)
        seq = seq[ind]
        
        dat = dat[ind,]
        
        
        # ----- remove too short peptides (l<6aa) -----
        l = nchar(seq)
        ind = which(l<3)
        if(length(ind)>0){dat = dat[-ind,]}
        
        dat$pep_seq = gsub("I","L",as.vector(dat$pep_seq))
        
        
        # ----- check which sequences could also be PCP and reassign them -----
        s = gsub("I","L",as.vector(dat$pep_seq))
        a = as.vector(dat$prot_acc)
        
        reassign <- numeric()
        
        for(i in 1:length(s)){
            
            k = grep("PSP_",a[i])
            if(length(k)>0){
                k2 = grep(s[i],proteome)
                if(length(k2)>0){reassign = c(reassign,i)}
            }
            
        }
        
        if(length(reassign)>0){
            
            for(i in 1:length(reassign)){
                
                k = strsplit(proteome,split=s[reassign[i]])[[1]][1]
                start = nchar(k)+1
                end = start+nchar(s[reassign[i]])-1
                a[reassign[i]] = paste("PCP_",start,"_",end,sep="")
                
            }
        }
        
        dat$prot_acc = a
        
        
        # ----- keep onyly highest score query -----
        
        rem = numeric()
        ambiguous <- list()
        ambiguous2 <- list()
        counter = 0
        counter2 = 0
        
        
        for(i in 1:max(dat$pep_query)){
            
            ind = which(dat$pep_query==i)
            
            if(length(ind)>1){

                d = dat[ind,]
                ind2 = which(d$pep_score==max(d$pep_score))
                
                if(length(ind2)==1){
                    
                    k = which(grepl(PSPstring,d$prot_acc[ind2])==TRUE)
                    
                    if(length(k)==0){
                        rem = c(rem,ind[-ind2])
                    }
                    if(length(k)>0){
                        
                        j = which(grepl(PCPstring,d$prot_acc)==TRUE)
                        if(length(j)>0){
                            
                            # compute difference in PSP top score and best PCP score
                            sPCP = max(d$pep_score[j])
                            sPSP = d$pep_score[ind2]
                            diff = 1-sPCP/sPSP
                            if(diff<thrDifpcp){
                                # keep top PCP
                                jj = which(d$pep_score[j]==sPCP)[1]
                                rem = c(rem,ind[-j[jj]])
                            }
                            else{
                                # compute difference in PSP top score and next best PSP score
                                sPSP2 = max(d$pep_score[-ind2])
                                sPSP = d$pep_score[ind2]
                                diff = 1-sPSP2/sPSP
                                if(diff<thrDif){
                                    rem = c(rem,ind)
                                }
                                else{
                                    rem = c(rem,ind[-ind2])
                                }
                                
                            }
                            
                        }
                        else{
                            # compute difference in PSP top score and next best PSP score
                            sPSP2 = max(d$pep_score[-ind2])
                            sPSP = d$pep_score[ind2]
                            diff = 1-sPSP2/sPSP
                            if(diff<thrDif){
                                rem = c(rem,ind)
                            }
                            else{
                                rem = c(rem,ind[-ind2])
                            }
                        }
                    }
                    
                }
                
                
                if(length(ind2)>1){
                    
                    k = which(grepl(PCPstring,d$prot_acc[ind2])==TRUE)
                    #print(PCPstring2)
                    if(length(k)==1){
                        k = k[1]
                        ind2 = ind2[k]
                        rem = c(rem,ind[-ind2])
                    }
                    
                    if(length(k)>1){
                        #   print("YES")
                        #print(PCPstring2)
                        k2 = which(grepl(PCPstring,d$prot_acc[ind2[k]])==TRUE)
                        if(length(k2)>0){
                            ind2 = ind2[k2[1]]
                            rem = c(rem,ind[-ind2])
                        }
                        else{
                            s = as.vector(d$pep_seq[ind2[k]])
                            # replace I and L by X
                            sm1 = gsub("I","X",s)
                            sm2 = gsub("L","X",sm1)
                            su = unique(sm2)
                            
                            if(length(su)==1){
                                k = k[1]
                                ind2 = ind2[k]
                                rem = c(rem,ind[-ind2])
                            }
                            if(length(su)>1){
                                counter2 = counter2+1
                                ambiguous2[[counter2]] = d[ind2,]
                                rem = c(rem,ind)
                            }
                        }
                        
                    }
                    if(length(k)==0){
                        
                        s = as.vector(d$pep_seq[ind2])
                        # replace I and L by X
                        sm1 = gsub("I","X",s)
                        sm2 = gsub("L","X",sm1)
                        su = unique(sm2)
                        
                        if(length(su)==1){
                            # ind2 = ind2[1]
                            # rem = c(rem,ind[-ind2])
                            
                            j = which(grepl(PCPstring,d$prot_acc)==TRUE)
                            if(length(j)>0){
                                
                                # compute difference in PSP top score and best PCP score
                                sPCP = max(d$pep_score[j])
                                sPSP = d$pep_score[ind2[2]]
                                diff = 1-sPCP/sPSP
                                if(diff<thrDifpcp){
                                    # keep top PCP
                                    jj = which(d$pep_score[j]==sPCP)[1]
                                    rem = c(rem,ind[-j[jj]])
                                    rem = c(rem,ind)
                                    
                                }
                                else{
                                    # compute difference in PSP top score and next best PSP score
                                    sPSP2 = max(d$pep_score[-ind2])
                                    sPSP = d$pep_score[ind2[2]]
                                    diff = 1-sPSP2/sPSP
                                    if(diff<thrDif){
                                        rem = c(rem,ind)
                                    }
                                    else{
                                        rem = c(rem,ind[-ind2])
                                    }
                                    
                                }
                                
                            }
                            else{
                                # compute difference in PSP top score and next best PSP score
                                sPSP2 = max(d$pep_score[-ind2])
                                sPSP = d$pep_score[ind2[1]]
                                diff = 1-sPSP2/sPSP
                                if(diff<thrDif){
                                    rem = c(rem,ind)
                                }
                                else{
                                    rem = c(rem,ind[-ind2[1]])
                                }
                            }
                            
                            
                        }
                        if(length(su)>1){
                            counter = counter+1
                            ambiguous[[counter]] = d[ind2,]
                            rem = c(rem,ind)
                        }
                    }
                }
            }
        }
        
        
        if(length(rem)>0){
            dat = dat[-rem,]
        }
        print("DONE")
        
        
        #----- filter for ion score and q-value -----
        
        cutOff = KK
        
        ind = which(dat$pep_score>IONscore)
        dat = dat[ind,]
        
        ind = which(dat$pep_expect>KK)
        if(length(ind)>0){dat = dat[-ind,]}
        
        
        
        #----- assign PCPs and PSPs -----
        # PCP
        ind = which(grepl(PCPstring,dat$prot_acc)==TRUE)
        pcp = dat[ind,]
        
        # PSPall
        ind = which(grepl(PSPstring,dat$prot_acc)==TRUE)
        psp = dat[ind,]
        
        #----- summarising peptides and types -----
        
        # counting peptides and printing summary
        s_psp = unique(as.vector(psp$pep_seq))
        s_pcp = unique(as.vector(pcp$pep_seq))
        
        
        print(paste0("PCP total ",length(s_pcp)," (",round(100*length(s_pcp)/(length(s_psp)+length(s_pcp)),2)," %)"))
        print(paste0("PSP total ",length(s_psp)," (",round(100*length(s_psp)/(length(s_psp)+length(s_pcp)),2)," %)"))
        print("............")
        
        #----- store extracted sequences -----
        save(pcp,file=paste0("OUTPUT/tmp/",cnt$protein_name[c],"/",indir,"/pcp.RData"))
        save(psp,file=paste0("OUTPUT/tmp/",cnt$protein_name[c],"/",indir,"/psp.RData"))
        
    }
    
    ### OUTPUT ###
    p = cnt$protein_name[cnt$substrateSeq == si][1]
    write("peptides were extracted",
          file = paste0("OUTPUT/tmp/",p,"/extract_peptides.txt"))
    
}


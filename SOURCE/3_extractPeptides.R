### qiSPI ###
# description:  extract peptides from search result file
#               - format MSDB according to invitroSPI
#               - re-map peptides to account for I/L redundancy
#               - apply PSM filtering / rescoring algorithm of invitroSPI
# input:        sample list (protein sequences), parsed search results
# output:       PCPs and PSPs per replicate
# author:       JL, HPR

library(seqinr)
library(stringr)
library(dplyr)

source("SOURCE/qiSPI_utils.R")

print("--------------------------------------------")
print("3) EXTRACT PEPTIDES FROM SEARCH RESULT FILES")
print("--------------------------------------------")

protein_name = snakemake@params[["protein_name"]]
PROTEINSEQUENCES_PATH = "INPUT/sequences/"

q_value = snakemake@params[["q_value"]]
ion_score = snakemake@params[["ion_score"]]
delta_score = snakemake@params[["delta_score"]]

lenCutoff = 5  # minimal product length


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
        
        # ----- format search results acc to invitroSPI format -----
        dat = dat %>%
            dplyr::rename(charge = pep_exp_z,
                          rank = pep_rank,
                          pepSeq = pep_seq,
                          ionScore = pep_score,
                          qValue = pep_expect,
                          PTM = pep_var_mod,
                          scanNum = pep_query)
        
        
        
        # ----- re-assign product types considering I/L redundancy -----
        #extract product type and position information
        pepName <- as.character(dat$prot_acc)
        pepType <- c()
        pepPos <- c()
        for (z in 1:length(pepName)){
            splitPepName <- unlist(strsplit(pepName[z], "_"))
            pepType <- c(pepType, splitPepName[1])
            splitPepName <- splitPepName[-1]
            pepPos <- c(pepPos, paste(splitPepName, collapse = "_"))
            
        }
        
        dat$productType <- pepType
        dat$spliceType = NA
        dat$positions <- pepPos
        
        # add metainfo
        dat$substrateID = cnt$substrateID[c]
        dat$substrateSeq <- proteome
        dat$digestTime <- cnt$digestTime[c]
        dat$runID = gsub(".raw","",cnt$raw_file[c])
        
        # assign substrate hits as PCP
        substratehits = which(dat$substrateSeq == dat$pepSeq)
        if (length(substratehits) > 0) {
            dat$productType[substratehits] = "PCP"
            dat$positions[substratehits] = paste(1, nchar(dat$substrateSeq[substratehits]), sep = "_")
        }
        
        print("")
        paste0("obtain correct annotations for positions and splice types for run: ",
               dat$runID[1]) %>%
            print()
        
        dat = dat %>% mapping()
        dat$spliceType[dat$productType == "PCP"] = NA
        
        # ----- remove empty queries and too short peptides -----
        seq = as.vector(dat$pep_seq)
        ind = which(is.na(seq) | nchar(seq)<lenCutoff)
        if(length(ind)>0){dat = dat[-ind,]}
        
        
        #----- filtering and re-scoring of PSMs -----
        
        PSMs = PSMfiltering(dat, delta_score = delta_score,
                            ion_score = ion_score, q_value = q_value)
        allPSMs = PSMs$allPSMs
        
        #----- assign PCPs and PSPs -----
        # PCP
        pcpidx = which(allPSMs$productType == "PCP")
        pspidx = which(allPSMs$productType == "PSP")
        
        pcp = allPSMs[pcpidx,]
        psp = allPSMs[pspidx,]
        
        #----- summarising peptides and types -----
        
        # counting peptides and printing summary
        s_psp = unique(as.vector(psp$pepSeq))
        s_pcp = unique(as.vector(pcp$pepSeq))
        
        
        print(paste0("unique PCPs: ",length(s_pcp)," (",round(100*length(s_pcp)/(length(s_psp)+length(s_pcp)),2)," %)"))
        print(paste0("unique PSPs: ",length(s_psp)," (",round(100*length(s_psp)/(length(s_psp)+length(s_pcp)),2)," %)"))
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


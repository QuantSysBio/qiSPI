### iqSPI ###
# description:  split quantitation results according to replicates
# input:        sample list, table matches
# output:       quantitation results for each time points
# author:       JL, modified by HPR

library(dplyr)
library(XML)
library(gtools)
library(stringr)

print("-----------------------------------------------------")
print("2) SPLIT QUANTITATION RESULTS ACCORDING TO REPLICATES")
print("-----------------------------------------------------")

protein_name = snakemake@params[["protein_name"]]
QUANTITATIONRESULTS_PATH = "INPUT/quantitation_results/"


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name %in% protein_name, ]

load(snakemake@input[["rawNames"]])


### MAIN PART ###
qresults = sample_list$quantitation_result %>%
    unique()

for (q in qresults) {
    
    print("")
    paste0("LOADING: ", q) %>% print()
    
    qfile = paste0(QUANTITATIONRESULTS_PATH, q)
    if (grepl(".txt",q)) {
        d = read.table(qfile,header=TRUE)
    } else if (grepl(".csv",q)) {
        d = read.csv(qfile, stringsAsFactors = F)
    } else if (grepl(".html",q)) {
        d = readHTMLTable(qfile, header = TRUE, as.data.frame = TRUE, which = 2)
    }
    
    # remove NA lines if any
    k = which(d$Component=="")
    if(length(k)>0){ d = d[-k,] }
    
    # sort file indices
    FileIndex = as.vector(d$Component)
    # FileIndexUniqueSorted = FileIndex %>%
    #     unique() %>%
    #     mixedsort()
    # FileIndexUniqueSorted = c("Ref", FileIndexUniqueSorted[grepl("C[:digit:]*", FileIndexUniqueSorted)])
    
    # split files and store in correct folder
    cnt = sample_list[sample_list$quantitation_result == q,]
    
    pb = txtProgressBar(min = 0, max = nrow(cnt), style = 3)
    for(i in 1:nrow(cnt)){
        setTxtProgressBar(pb, i)
        
        k = which(FileIndex == cnt$order[i])
        dat = d[k,]
        
        cntRawNames = rawNames[names(rawNames) == cnt$protein_name[1]] %>%
            unlist()
        names(cntRawNames) = str_split(names(cntRawNames), coll("."), simplify = T)[,2] %>%
            as.character()
        outdir = names(cntRawNames)[cntRawNames == cnt$order[i] & names(cntRawNames) %in% gsub(".raw","",cnt$raw_file)]
        write.csv(x=dat,
                  file=paste0("OUTPUT/tmp/",cnt$protein_name[i],"/",outdir,"/quantitationResults.csv"),
                  row.names = F)
        
    }
    
    ### OUTPUT ###
    p = sample_list$protein_name[sample_list$quantitation_result == q][1]
    write("quantitation results were splitted",
          file = paste0("OUTPUT/tmp/",p,"/split_quantitation.txt"))
    
}





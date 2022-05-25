### iqSPI ###
# description:  - remove headers from search result files
#               - create temporary output directories
# input:        sample list
# output:       folder structure with search results for each replicate
# author:       JL, modified by HPR

library(dplyr)

print("--------------------------------------------------------")
print("1) PARSE SEARCH RESULT FILES AND CREATE FOLDER STRUCTURE")
print("--------------------------------------------------------")

protein_name = snakemake@params[["protein_name"]]
SEARCHRESULTS_PATH = "INPUT/search_results/"


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       stringsAsFactors = F)
# sample_list = read.csv(file = "INPUT/sample_list.csv",
#                        stringsAsFactors = F)

sample_list = sample_list[sample_list$protein_name %in% protein_name, ]

### MAIN PART ###
if (! paste0("OUTPUT/tmp/", protein_name) %>% dir.exists() %>% all()) {
  sapply(paste0("OUTPUT/tmp/", protein_name), dir.create)
}

numQueries = list()
rawNames = list()

prots = sample_list$protein_name %>%
  unique()

for (pp in 1:length(prots)) {
  
  p = prots[pp]
  print("")
  paste0("LOADING: ", p) %>% print()
  
  cnt = sample_list[sample_list$protein_name == p,]
  cntNumQueries = rep(NA, nrow(cnt))
  cntRawNames = rep(NA, nrow(cnt))
  
  pb = txtProgressBar(min = 0, max = nrow(cnt), style = 3)
  for (i in 1:nrow(cnt)) {
    setTxtProgressBar(pb, i)
    
    search_results = read.csv(paste0(SEARCHRESULTS_PATH, cnt$search_result[i]),
                              stringsAsFactors = F)
    protein = cnt$protein_name[i]
    raw = cnt$raw_file[i]
    
    # create temporary folder for search results per replicate / time point
    outdir = paste0("OUTPUT/tmp/", protein, "/", gsub(".raw","",raw))
    if (!dir.exists(outdir)) {
      dir.create(outdir)
    }
    
    # load search result file
    con = file(paste0(SEARCHRESULTS_PATH, cnt$search_result[i]),"r")
    i1 = readLines(con)
    close(con)
    
    # remove header from search result file
    start = grep("prot_hit_num",i1)
    end = grep("Peptide matches not assigned to protein hits",i1)-1
    
    tmp = i1[grep("Number of queries",i1)]
    cntNumQueries[i] = as.numeric(strsplit(tmp,split=",")[[1]][2])
    names(cntNumQueries)[i] = gsub(".raw","",raw)
    
    cntRawNames[i] = cnt$order[i]
    names(cntRawNames)[i] = gsub(".raw","",raw)
    
    if(length(end) > 0) {
      search_results = i1[start:end]
    } else {
      search_results = i1[start:length(i1)]
    }
    
    write(search_results,
          file = paste0(outdir, "/searchResults.csv"),
          ncolumns=1)
    
  }
  
  numQueries[[pp]] = cntNumQueries
  rawNames[[pp]] = cntRawNames
  names(numQueries)[pp] = protein
  names(rawNames)[pp] = protein
  
}

### OUTPUT ###

save(numQueries, file = unlist(snakemake@output[["numQueries"]]))
save(rawNames, file = unlist(snakemake@output[["rawNames"]]))


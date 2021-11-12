### iqSPI ###
# description:  split quantities according to biological replicates
# input:        sample list, merged quantification results
# output:       quantities (merged over replicates)
# author:       HPR, adapted from JL


library(dplyr)
library(stringr)

print("-------------------------------------------------------")
print("7) SPLIT INTENSITIES ACCORDING TO BIOLOGICAL REPLICATES")
print("-------------------------------------------------------")

protein_name = snakemake@params[["protein_name"]]


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name %in% protein_name, ]


### MAIN PART ###

prots = sample_list$protein_name %>%
  unique()

for (p in prots) {
  
  paste0("merging replicates of ", p) %>%
    print()
  
  # ----- load merged intensities -----
  load(paste0("OUTPUT/",p,"/quantity.RData"))
  
  # ----- sort by sequence -----
  k = order(quantity[,2])
  quantity = quantity[k,]
  
  c = which(names(quantity) == "z")
  
  # ----- split replicates -----
  cnt = sample_list[sample_list$protein_name == p,]
  reps = cnt$biological_replicate %>%
    unique() %>%
    sort(decreasing = F)
  
  rep = list()
  for (r in 1:length(reps)) {
    
    x = which(cnt$biological_replicate == reps[r])
    q = quantity[, c(1:c, c+x)]
    names(q) = c("protHit", "sequence", "z", paste0("tp_",
                                                    str_replace_all(as.character(cnt$digestTime[x]),coll("."),coll(","))))
    rep[[r]] = q
    
  }
  
  names(rep) = reps
  
  # ----- save replicates -----
  save(rep, file=paste0("OUTPUT/",p,"/quantity_rep.RData"))
  write.csv(rep, file=paste0("OUTPUT/",p,"/quantity_rep.csv"),
            row.names = F)
  
}


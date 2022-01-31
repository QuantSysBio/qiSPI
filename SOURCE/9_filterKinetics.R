## qiSPI ###
# description:  plot processed kinetics, aggregate over mean per replicate
# input:        results.RData
# output:       plotted kinetics
# author:       JL, modified by HPR


library(dplyr)
library(RColorBrewer)

print("-------------------------------------")
print("9) FILTERING AND PLOTTING OF KINETICS")
print("-------------------------------------")

protein_name = snakemake@params[["protein_name"]]


### INPUT ###
sample_list = read.csv(file = snakemake@input[["sample_list"]],
                       stringsAsFactors = F)
sample_list = sample_list[sample_list$protein_name %in% protein_name, ]


### MAIN PART ###
# ----- colour definition -----
c1 <- colorRampPalette(brewer.pal(8, "Set1"))(9)
c2 <- colorRampPalette(brewer.pal(8, "Set2"))(8)
c3 <- colorRampPalette(brewer.pal(8, "Set3"))(12)
c4 <- colorRampPalette(brewer.pal(8, "Paired"))(12)
c5 <- colorRampPalette(brewer.pal(8, "Accent"))(8)

myColors = c4[c(1,3,4,6,7,9)]

# --- load results, filter and plot -----

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
    
    
    # ----- load results.RData and plot unfiltered -----
    
    load(file=paste0("OUTPUT/",p,"/results.RData"))
    outdir = paste0("OUTPUT/",p,"/PLOTS/")
    
    if (!dir.exists(outdir)) {
        dir.create(outdir)
    }
    
    pdf(paste0(outdir,"unfiltered.pdf"), width=10, height=10)
    par(mfrow = c(4,5))
    
    for(i in 1:dim(results[[1]])[1]){
        
        maxi = 0
        mini = 10**20
        for(j in 1:length(results)){
            maxi = max(c(maxi,as.numeric(as.vector(results[[j]][i,2:(1+length(t))]))))
            mini = min(c(mini,as.numeric(as.vector(results[[j]][i,2:(1+length(t))]))))
        }
        
        plot(1,1,type="b",col="white",xlim=c(0,max(t)),axes=FALSE,main=paste(i,results[[j]][i,dim(results[[j]])[2]],sep="_"),ylim=c(mini,maxi),xlab="time",ylab="Intensity")
        axis(1)
        axis(2)
        
        for(j in 1:length(results)){
            points(t,results[[j]][i,2:(1+length(t))],type="b",col=myColors[j])
        }
        
    }
    
    dev.off()
    
    
    # ----- subtract noise -----
    for (r in 1:length(results)) {
        df = results[[r]][, c(2:(1+length(t)))]
        
        # subtract minimum intensity
        tmp = apply(df, 1, function(x) {
            x = as.numeric(x)
            k = which.min(x)
            
            x = x - x[k]
            return(x)
            
        }) %>%
            t() %>%
            as.data.frame()
        
        names(tmp) = names(df) 
        
        # set 0 hour tp to 0
        z = which(names(tmp) %in% c(0, "CTRL"))
        if (length(z) > 0) {
            tmp[,z] = 0
        }
        
        results[[r]][, c(2:(1+length(t)))] = tmp
    }
    
    # ----- get mean over technical replicates -----
    bioReps = cnt$biological_replicate[cnt$digestTime == t[1]] %>% 
        unique()
    
    M = list()
    for (b in bioReps) {
        M[[b]] = matrix(NA,dim(results[[1]])[1],length(t))
    }
    
    pdf(paste0(outdir,"unfiltered_means.pdf"), width=10, height=10)
    par(mfrow = c(4,5))
    
    for(i in 1:dim(results[[1]])[1]){
        
        maxi = 0
        mini = 10**20
        for(j in 1:length(results)){
            maxi = max(c(maxi,as.numeric(as.vector(results[[j]][i,2:(1+length(t))]))))
            mini = min(c(mini,as.numeric(as.vector(results[[j]][i,2:(1+length(t))]))))
        }
        
        plot(1,1,type="b",col="white",xlim=c(0,max(t)),
             axes=FALSE,
             main=paste(i,results[[j]][i,dim(results[[j]])[2]],sep="_"),
             ylim=c(mini,maxi),xlab="time",ylab="Intensity")
        
        axis(1)
        axis(2)
        
        for (b in bioReps) {
            
            bb = which(cnt$biological_replicate[cnt$digestTime == t[1]] == b)
            x = sapply(results[bb], function(arg) {
                as.numeric(as.vector(arg[i,2:(1+length(t))]))
            }) %>%
                as.matrix() %>%
                t()
            
            m = apply(x,2,mean)
            points(t,m,type="b",col=myColors[2*b-1])
            M[[b]][i,] = m
            
        }
    }
    
    
    dev.off()
    
    # ----- sort by max intensity detected -----
    maxi = rep(NA,dim(M[[1]])[1])
    for(i in 1:dim(M[[1]])[1]){
        
        x = numeric()
        for(j in 1:length(M)){
            x = c(x,M[[j]][i,])
        }
        maxi[i] = max(x)
        
    }
    
    ind = order(maxi,decreasing = TRUE)
    for(j in 1:length(M)){
        results[[j]] = results[[j]][ind,]
        M[[j]] = M[[j]][ind,]
    }
    
    # ----- save filtered results -----
    
    filteredResults = results
    filteredMeans = M
    save(filteredResults,
         file=paste0("OUTPUT/",p,"/filteredResults.RData"))
    save(filteredMeans,
         file=paste0("OUTPUT/",p,"/filteredMeans.RData"))
    
    
}


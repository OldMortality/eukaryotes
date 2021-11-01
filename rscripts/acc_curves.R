rm(list=ls())
setwd("C:/Users/michel/Documents/eukaryotes/rscripts")
library(dplyr)
library(edgeR)


source("euk_functions.R")

set.seed(133499)


makeDF <- function (d,sample) {
  dd <- d[which(d$Sample==sample),c("Abundance","phylum")]
  d3 <- dd %>% 
    group_by(phylum) %>% 
    summarise(occs = sum(Abundance)) %>% as.data.frame()
  n <- sum(d3$occs)
  r <- character()
  for (i in 1:nrow(d3)) {
    z <- rep(d3[i,'phylum'],d3[i,'occs'])
      r <- c(r,z)
  }
  stopifnot(length(r)==n)
  return(r)
}

#read_depth <- 500



makeCurveData <- function(sampleDF) {
  curve_data <- vector()
  s <- sample(1:length(sampleDF),replace=F)
  reads <- sampleDF[s]
  
  for (read_depth in seq(1,length(sampleDF))) {
    n <- length(unique(reads[1:read_depth]))
    curve_data[read_depth] <- n
    if (n == length(unique(sampleDF))) {
      break
    }
  }
  return(curve_data)
}

calcTuring <- function(d, sample) {
  dd <- d[which(d$Sample==sample),c("Abundance","phylum")]
  d3 <- dd %>% 
    group_by(phylum) %>% 
    summarise(occs = sum(Abundance)) %>% as.data.frame()
  freqs <- d3$occs
  t <- goodTuring(freqs)
  return(t$P0)
}




makePlots <- function(sampleDF, sample, N_SIMS, counter) {
  

  c1 <- makeCurveData(sampleDF)

  p1 <- ggplot(data=data.frame(x=seq(1,length(c1)),y=c1)) + 
                              geom_point(aes(x=x,y=y)) + 
    xlab("Number of OTUs read") + 
    ylab("Number of phyla discovered") +
    ggtitle(paste0(counter,'. ',sample))
  #p1
  
  # mininum number of reach for each simulation to find all phyla
  min_sims <- numeric()
  
  for (i in 1:N_SIMS) {
    c1 <- makeCurveData(sampleDF)
    
    # the number of reads which first found all phyla
    min_reads <- min(which(c1==length(unique(sampleDF))))
    min_sims[i] <- min_reads
  }
  
  p2 <- ggplot(data=data.frame(x=min_sims)) + geom_histogram(aes(x=x),colour='black', fill='white') +
    xlab("# OTUs needed to find all phyla in the sample") +
    ggtitle(sample)
  
  highest_sim <- max(min_sims)
  p95 <- as.numeric(quantile(min_sims,0.95))
  tu <- calcTuring(d = d,sample = sample)  
  
  return(list(p1=p1,p2=p2,n_phyla=length(unique(sampleDF)),highest_sim=highest_sim,p95=p95,tu=tu))
}



d <- loadData()
samples <- sort(unique(d$Sample))
#samples <- unique(d$Sample)

phyla <- sort(unique(d$phylum))
d <- d[-which(d$phylum=="no hit, reference, or complete taxonomy string"),]

# sample <- samples[1]
# sampleDF <- makeDF(d,sample)

counter <- 1
results <- data.frame()
pdf('report.pdf',width = 14,height=8)
for (sample in samples) {
  sampleDF <- makeDF(d,sample)
  r <- makePlots(sampleDF = sampleDF, sample=sample, N_SIMS = 1000,counter=counter)
  p1 <- r$p1
  p2 <- r$p2
  
  g <- grid.arrange(p1,p2,nrow=1)
  new_row <- c(sample,r[["n_phyla"]],r[["highest_sim"]],r[["p95"]],r[["tu"]])
  results <- rbind(results,new_row)
  print(Sys.time())
  print(new_row)
  counter <- counter + 1
  
}
dev.off()
colnames(results) <- c('Sample',"n phyla","Highest_Sim","p95","P0")

get_total_abundance_by_sample <- function(df) {
  dd <- df[,c("Sample","Abundance","phylum")]
  result <- dd %>% 
    group_by(c(Sample)) %>% 
    summarise(occs = sum(Abundance)) %>% as.data.frame()
  colnames(result) <- c("Sample","Sample_Abundance")
  return(result)
}

total_abundance_by_sample <- get_total_abundance_by_sample(df=d)
results$sample_abundance<- total_abundance_by_sample[match(results$Sample,total_abundance_by_sample$Sample),'Sample_Abundance']
results$p95 <- ceiling(as.numeric(results$p95))
results$`n phyla` <- as.numeric(results$`n phyla`)
results$sample_abundance <- as.numeric(results$sample_abundance)
results$Highest_Sim <- as.numeric(results$Highest_Sim)
results
save.image(file="after_loop.rdata")

results

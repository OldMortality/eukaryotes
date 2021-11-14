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
results <- results[-which(results$`n phyla`==1),]
results$frac <- round(results$p95 / results$sample_abundance,2)

#save.image(file="after_loop.rdata")
#load(file="after_loop.rdata")


####
#### Absorbing Markov Chain Approach
####

library(binaryLogic)


new_ph <- function(x,y,p,N) {
  result <- 0
  #print(paste(x,y))
  if (x==y) return(NA) 
  b1 <- as.numeric(as.character(as.binary(x-1,n = N)))
  b2 <- as.numeric(as.character(as.binary(y-1,n = N)))
  if ( ( sum(b2) == sum(b1) + 1 ) && (length(which(b1!=b2)) == 1)  && (b2[which(b1!=b2)]==1) ) {
    result <- p[which(b1!=b2)]
  }
  return(result)
}




calcStats <- function(sampleId, p) {
  N <- length(p) # number of phyla
  T <- matrix(0,nrow=2^N,ncol=2^N)
  for (i in seq(1,2^N-1)) {
    T[i,i] <- sum(p * as.numeric(as.character(as.binary(i-1,n = N)))) 
  }
  T[2^N,2^N] <- 1
  
  for (i in seq(1,2^N-1)) {
    for (j in seq(i + 1,2^N)) {
      if (i != j) {
        T[j,i] <- new_ph(x=i,y=j,p=p,N=N)  
      }
      
    }
  }
  T
  
  Q <- t(T[1:2^N-1,1:2^N-1])
  I <- matrix(0,nrow=2^N-1,ncol=2^N-1)
  for (i in 1:2^N-1) {
    I[i,i] <- 1
  }
  
  M <- solve(I-Q)
  one_vector <- matrix(1,nrow=2^N-1,ncol=1) 
  # expected number of steps
  mn <- M %*% one_vector
  # variance of the expected number of steps
  s <- (2 *  M - I ) %*% mn - mn * mn
  ## initial state is 1
  mn <- mn[1,1]
  s <- s[1,1]
  return(list(mn = mn, s=s, upp = mn + 2 * sqrt(s)))
    
}

p <- c(0.2,0.1,0.7)

makeProps <- function(d,sample) {
  dd <- d[which(d$Sample==sample),c("Abundance","phylum")]
  d3 <- dd %>% 
    group_by(phylum) %>% 
    summarise(occs = sum(Abundance)) %>% as.data.frame()
  n <- sum(d3$occs)
  p <- d3$occs/n
  lowest <- p[which(p == min(p))]
  p <- sort(c(p,lowest/2))
  p <- p/sum(p)
  return(p)
}

dff <- NULL
for (s in 1:length(samples)) {
  print(s)
  c <- calcStats(sampleId = samples[s],p=makeProps(d=d,sample=samples[s]))
  dff <- rbind(dff,data.frame(sample=samples[s],mn=c[["mn"]],s=c[["s"]],upp=c[["upp"]]))
}
#save.image(file='after_mc.rdata')
load(file='after_mc.rdata')
m <- merge(results,dff,by.x='Sample',by.y='sample')
head(m)


missed <- m[which(m$mn>m$sample_abundance),]
nrow(missed)
dd <- d[which(d$Sample==missed$Sample[1]),c("Abundance","phylum")]
d3 <- dd %>% 
  group_by(phylum) %>% 
  summarise(occs = sum(Abundance)) %>% as.data.frame()
d3

#save.image(file='after_mc3.rdata')
save(m,file='m.rdata')


dd <- d[which(d$Sample==samples[1]),c("Abundance","phylum")]
d3 <- dd %>% 
  group_by(phylum) %>% 
  summarise(occs = sum(Abundance)) %>% as.data.frame()
d3

# tf <- data.frame(a=1:10,b=1:10,c=1:10)
# library(purrr)
# pmap_dfr(tf, function(a, b, c) {
#   data.frame(var1 = a + b * c,
#              var2 = c/2) 
# }  )

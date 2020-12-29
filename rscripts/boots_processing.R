##
## makes the plots from the bootstrap results.
##
rm(list=ls())
library(gridExtra)

setwd("C:/Users/Michel/Documents/eukaryotes/rscripts")
source('boots_processing_1.R')
load(file='varnames.RData')



phyls <- c('Arthropoda',
            'Bacillariophyta',
            'Cercozoa',
            'Chytridiomycota',
            'Nematoda',
            'Streptophyta',
           'Ascomycota',
           'Basidiomycota',
           'Chlorophyta',
           'Ciliophora',
           'Rotifera',
           'Tardigrada')

ps <- list()
for (i in 1:length(phyls)) {
  ph <- phyls[i]
  fname <- paste0('../boots/',ph,'.RData')
  load(file=fname)
  ps[[i]] <- process_bootstraps(phylum=ph,
                                bm2=bm2,
                                var.names=var.names)
}           

for (i in 1:length(phyls)) {
  
  p1a <- ps[[i]]$plot.zeros
  p1b <- ps[[i]]$plot.ci
  p.combined <- grid.arrange(p1a,p1b,ncol=2)
  ggsave(p.combined, file=paste0("../images/boots2_",phyls[i],".bmp"), 
                         device="bmp",dpi=800,
         width=200,height=100,units='mm')
  
}

# load(file="allstats.RData")
# all.stats[[2]]

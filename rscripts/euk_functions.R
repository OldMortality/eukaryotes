## 
## R code to analyse Paul's eukaryotes data from Antarctica
##
## This file contains all helper functions. The main ones are
##    * create.df.species(filepath,locations)
##        this one creates a dataframe with all relevant samples
##    * getDFByPhylum(df,phylum,cols)
##        creates a dataframe for your phylum. Do typically you would do:
##          df.species <- create.df.species(filepath=..,locations=LOCATIONS)
##        and then for each phylum:
##          df.tardigrades <- getDFByPhylum(df.species,'Tardigrada',LR.COLUMNS)
##        Now you can use df.tardigrades for logistic regression.
##
##
library(dplyr)
library(pipeR)


## Global constants
##
FILEPATH = 'C:/Users/Michel/Documents/eukaryotes/data/200_all_data_long_export_filtered.Rdata'
 

# Global variable for our 3 locations
LOCATIONS <- c('Lake Terrasovoje','Mawson Escarpment','Mount Menzies')

# All columns of interest
MYCOLNAMES <- c("AMMN","NITR","PHOS","POTA","SULPH","CARB","COND","PH_CACL",             
                "RLU","QUARTZ","FELDSPAR","TITANITE","GARNETS","MICAS",
                "DOLOMITE","KAOLCHLOR","CALCITE","CHLORITE","SLOPE")

# Columns used in logistic regression. I left out many of the NA columns.
LR.COLUMNS <- c("present", "Location", "Abundance",
                "POTA" , "SULPH" , "COND" , "PH_CACL" ,  
                "RLU" , "QUARTZ" , "FELDSPAR" , "TITANITE" , "GARNETS" , 
                "MICAS" , "DOLOMITE" , "KAOLCHLOR" , "CALCITE","CHLORITE")

## Inverse logistic function
invlogit <- function(x) {
  return(1/(1+exp(-x)))
}

# returns vector of sample-id's of those samples with low total abundance
getSamplesWithLowAbundance <- function(df,threshold=1000) {
  abundances    <- aggregate(Abundance ~ Sample, df, sum)
  low.abundance <- abundances[which(abundances$Abundance < threshold),c("Sample")]
  return(low.abundance)  
}
 
removeSamplesWithLowAbundance <- function(df,threshold=1000) {
  lows <- getSamplesWithLowAbundance(df,threshold)
  dropm <- which(df$Sample %in% lows)
  result <- df
  if (length(dropm > 0)) {
    result <- df[-dropm,] 
  }
  return(result)
}


# read data from file, and return as data.frame
loadData <- function(filepath=FILEPATH) {
  load(filepath)
  return(psob_molten)
}

# keepOnlyEukaryotes <- function(df) {
#   return(subset(df,superkingdom %in% 'Eukaryota'))
# }

# keepOnlyMyLocations <- function(df,myLocations) {
#   return(subset(df,Location %in% myLocations ))
# }


## get all that aren't there :)
# getAbsences <- function(df) {
#   return(which(df$Abundance==0))
# }

# removeAbsences <- function(df) {
#   abs <- getAbsences(df)
#   if (length(abs)==0) {
#     result <- df
#   } else {
#     result <- df[-getAbsences(df=df),]  
#   }
#   return(result)
# }

# Remove rows where we have the same species for multiple OTUs. 
#   If we do this, anything that follows will only work for presence/absence, but not for counts. 
# removeDupSpecies <- function(df) {
#   result <- df[!duplicated(df[,c('Sample','species')]),]
#   return(result)
# }


library(dplyr) # for join
create.df.species <- function(filepath=FILEPATH) {
  d1 <- loadData(filepath = filepath)   
  return(d1)
}

takeLogoffactor <- function(x) {
  return(log(as.numeric(as.character(x))))
}

# take the log of the soil data. These have been read in
#   as factors, so we need to do as.character() first
# takeLogs <- function(df,colnames) {
#   df[,colnames] <- apply(df[,colnames],2,FUN= log) #takeLogoffactor)
#   # we don't want logs of PH
#   df$PH_H2O <- exp(df$PH_H2O)
#   df$PH_CACL <- exp(df$PH_CACL)
#   return(df)
# }



# returns 1 row for each sample, containing soil data
#    
getSoilData <- function(df,colnames=MYCOLNAMES) {
  # Get the first of each Sample, and for those get the columns. They
  #   should be the same for each for the same sample
  s  <- df[!duplicated(df[,c('Sample')]),c(c("Sample","Location"),colnames)]
  return(s)
}

# get a dataframe, one row for each sample, and a column 'present',
#  which tells us whether the phylum is present in each sample
#  df is by species, so there could be many  rows in df with the phylum,
#  but only 1 row per sample is returned.
getSamplesWithPhylumPresence <- function(df,phylum) {
  df <- as.data.frame(df)
  all.samples <- unique(df$Sample)
  samples.with.phylum <- (unique(df[which(df$phylum==phylum),'Sample']))
  df.presence <- data.frame(Sample = all.samples)
  df.presence$present <- 0
  df.presence[which((df.presence$Sample) %in% samples.with.phylum),'present'] <- 1
  return(df.presence)
}


##
## df.species is a dataframe of all samples with positive abundance in any of our three locations. 
## If a species has multiple OTUs, it appears in this dataframe only once per sample.
##
# df.species <- create.df.species(filepath = '~/Documents/eukaryotes/data/200_all_data_long_export_filtered.Rdata')
# 

# return df with total abundance for each sample
#
getAbundancesBySample <- function(df) {
  df.abundances <- aggregate(Abundance ~ Sample, df, sum)
  return(df.abundances)
}

# for a given phylum, list how many distinct otu's there are in each sample
getNumberOfDistinctOTUSbySample <- function(df,phylum) {
  #& df$species=="Mesobiotus furciger"
  df2 <- df[which(df$phylum==phylum),]
  z <- tapply(df2$OTU, df2$Sample, FUN = function(x) length(unique(x)))
  df.z <- data.frame(Sample=rownames(z),distinct.otus=z)
  # merge in the samples with zero OTUs of this phylum
  all.samples <- data.frame(Sample=unique(df$Sample))
  m <- merge(all.samples,df.z,by='Sample',all.x=T)
  m[which(is.na(m$distinct.otus)),'distinct.otus'] <- 0
  return(m)
  
}

## Get all information we need for logistic regression for a given phylum
## Sample, Location, soildata, total abundance (in the sample), present/absent 
##   df will typically be df.species, and phylum would be e.g. 'Tardigrada'
getAllByPhylum <- function(df,phylum,colnames=MYCOLNAMES) {
  print(phylum)
  df.distinctOTUs <- getNumberOfDistinctOTUSbySample(df,phylum)[,c('Sample','distinct.otus')]   
  df.soilData  <-  getSoilData(df,colnames = colnames)
  result <- (merge(df.distinctOTUs,df.soilData,by='Sample'))
  result$SLOPE <- as.numeric(result$SLOPE)
  result$PH_CACL <- as.numeric(result$PH_CACL)
  numCols <- LR.COLUMNS[-c(1,2,3)]
  # Convert soildata from char to numeric
  for (i in 1:length(numCols)) {
    result[,numCols[i]] <- as.numeric(result[,numCols[i]]) 
  } 
  dim(result)
  return(result)
}


lineUp <- function(m,d) {
  par(mfrow=c(4,5))
  for (i in 1:19) {
    sim <- unlist(simulate(m))
    leaveOut <- as.numeric(attributes(m$na.action)$names)
    if (length(leaveOut > 0)) {
      d <- d[-leaveOut,]  
    }
    
    m2 <- glm(sim ~  log(Abundance) +  FELDSPAR +  MICAS ,data=d,family='binomial')
    hist(residuals(m2),probability = T)
  }
  hist(residuals(m),probability = T)
  
} 


scaleV <- function(v){
  ma <- max(v,na.rm=T)
  mi <- min(v,na.rm=T)
  v <- (v - mi)/(ma-mi)
  return(v)
}

makedummy <- function(vec, makenames = TRUE, contrasts = FALSE) {
  z <- unique(vec)
  X <- matrix(0, nrow = length(vec), ncol = length(z))
  X[cbind(1 : length(vec), match(vec, z))] <- 1
  if (makenames) colnames(X) <- paste0(deparse(substitute(vec)), "_", z)
  if (contrasts) X <- X[, -ncol(X)]
  return(X)
}





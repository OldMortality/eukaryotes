library(RUnit)

saveLoadedData <- NULL
ps_molten_saved <- NULL

test.invlogit <- function() {
  checkEquals(invlogit(0), 0.5)
}

test.loadData <- function() {
  path <- '~/Documents/eukaryotes/data/200_all_data_long_export_filtered.Rdata'
  df <- loadData(filepath=path)
  ps_molten_saved <- df
  checkEquals(dim(df)[1],1351680)
  checkEquals(dim(df)[2],48)
}

test.getSamplesWithLowAbundance <- function() {
  d <- data.frame(Sample=c(1,1,2,2,3),
                  Abundance=c(700,200,800,300,1200))
  
  lows <- getSamplesWithLowAbundance(d,threshold=1000)
  checkEquals(lows,c(1))
  ##
  d <- data.frame(Sample=c(1,1,2,2,3),
                  Abundance=c(700,200,800,300,950))
  
  lows <- getSamplesWithLowAbundance(d,threshold=1000)
  checkEquals(lows,c(1,3))
}

test.removeSamplesWithLowAbundance <- function() {
  d <- data.frame(Sample=c(1,1,2,2,3),
                  Abundance=c(700,200,800,300,950))
  d2 <- removeSamplesWithLowAbundance(d,threshold = 1000)
  correctAnswer <- d[which(d$Sample==2),]
  checkEquals(d2 ,correctAnswer)
  ##
  d[5,'Abundance'] <- 1100
  d2 <- removeSamplesWithLowAbundance(d,threshold = 1000)
  correctAnswer <- d[which(d$Sample==2 | d$Sample==3),]
  checkEquals(d2 ,correctAnswer)
}

test.keepOnlyEukaryotes <- function() {
  d <- data.frame(Sample=seq(11,20),
                  superkingdom=rep(c('Eukaryota','Bacteria'),5))
  d2 <- keepOnlyEukaryotes(d)
  correctAnswer <- d[which(d$superkingdom == 'Eukaryota'),]
  checkEquals(d2 ,correctAnswer)
}

test.getAbsences  <- function() {
  d <- data.frame(Sample=seq(1,100),
                  Abundance=round(rnorm(100,mean=1000)))
  zeroes <- c(10,20,30,40,50,100)
  d[zeroes,"Abundance"] <- 0
  checkEquals(zeroes,getAbsences(d))
}

test.removeAbsences <-  function() {
  d <- data.frame(Sample=seq(1,100),
                  Abundance=round(rnorm(100,mean=1000)))
  zeroes <- c(10,20,30,40,50,100)
  d[zeroes,"Abundance"] <- 0
  presences <- removeAbsences(d)
  checkEquals(length(which(presences$Abundance == 0)),0)
}

test.removeDupSpecies <- function() {
  d <- data.frame(Sample=rep('A',5),
                  species=c('tiger','tiger','lion','cheetah','leopard'))
  checkEquals(d[-2,],removeDupSpecies(d))
  d <- data.frame(Sample=rep('A',5),
                  species=c('tiger','tiger','lion','cheetah','tiger'))
  checkEquals(d[-c(2,5),],removeDupSpecies(d))
  d <- data.frame(Sample=rep('A',5),
                  species=c('tiger','tiger','lion','lion','tiger'))
  checkEquals(d[c(1,3),],removeDupSpecies(d))
  d <- data.frame(Sample=rep('A',5),
                  species=c('tiger','ocelot','lion','leopard','lynx'))
  checkEquals(d,removeDupSpecies(d))
}

test.removeNoPhylum <- function() {
  d <- data.frame(Sample=seq(1,5),
                  phylum=c('Chordata','Tardigrade','undefined','Tardigrade','undefined'))
  d2 <- removeNoPhylum(d) 
  checkEquals(d2$Sample,c(1,2,4)) 
}

test.create.df.species <- function() {
  df.species <- create.df.species(filepath = '~/Documents/eukaryotes/data/200_all_data_long_export_filtered.Rdata',
                                  locations = LOCATIONS)
  checkEquals(length(which(df.species$Abundance==0)),0)
  checkEquals(length(which(df.species$superkingdom != 'Eukaryota')),0)
  checkEquals(dim(df.species)[1] > 0, TRUE)
  checkEquals(dim(df.species)[2] > 0, TRUE)
  # no duplicate species in sample
  checkEquals(sum(duplicated(df.species[,c('Sample','species')])),0)
  # nothing in the wrong locations
  checkEquals(length(which(df.species$Location %in% LOCATIONS)),dim(df.species)[1])  
  saveLoadedData <- df.species
  
}

test.takeLogoffactor <- function() {
  x = c("1","2.71827")
  checkEqualsNumeric(c(0,1),takeLogoffactor(x),tolerance=0.001)
}

test.takelogs <- function() {
  d <- data.frame(PH_H2O = c("1","1"),PH_CACL=c("1","1"),x=c("2.71827","10"))
  t <- data.frame(PH_H2O = c(1,1),PH_CACL=c(1,1),x=c(1,2.30))
  d2 <- takeLogs(d,colnames = c("PH_H2O","PH_CACL","x"))
  checkEqualsNumeric(d2,t,tolerance=0.001)  
}

getGetSoilData <- function() {
  df <- data.frame(Sample=c('A','A','B','B'),
                   Location = c('X','X','Y','Y'),
                   x = c('1','1','2.718282','2.718282'),
                   y = c('25','25','20','20'),
                   PH_CACL = c('1','1','10','10'),
                   PH_H2O = c('1','1','10','10'))
  soilData <- getSoilData(df,colnames=c('x','y','PH_CACL',"PH_H2O" ))
  target <- data.frame(Sample=c('A','B'),
                       Location=c('X','Y'),
                       x=c(0,1),
                       y=c(log(25),log(20)),
                       PH_CACL=c(1,10),
                       PH_H2O=c(1,10))
  checkEquals(soilData$Sample,target$Sample)
  checkEquals(soilData$Location,target$Location)
  checkEqualsNumeric(soilData$x,target$x,tolerance = 0.01)
  checkEqualsNumeric(soilData$y,target$y,tolerance = 0.01)
  
  checkEqualsNumeric(soilData$PH_CACL,target$PH_CACL,tolerance = 0.01)
  checkEqualsNumeric(soilData$PH_H2O,target$PH_H2O,tolerance = 0.01)
}

test.getGetSamplesWithPhylumPresence <- function() {
  df <- data.frame(Sample=c('A','B','C','A','B','C'),
                   phylum=c('a','b','a','b','b','d'))
  df.result <- getSamplesWithPhylumPresence(df,'a')
  target <- data.frame(Sample=c('A','B','C'),
                       present=c(1,0,1))
  checkEquals(target,df.result)
  ##
  df <- data.frame(Sample=c('A','B','C','A','B','C'),
                   phylum=c('a','c','a','b','b','d'))
  df.result <- getSamplesWithPhylumPresence(df,'b')
  target <- data.frame(Sample=c('A','B','C'),
                       present=c(1,1,0))
  checkEquals(target,df.result)
  ##
  df <- data.frame(Sample=c('A','B','C','A','B','C'),
                   phylum=c('a','e','a','b','e','d'))
  df.result <- getSamplesWithPhylumPresence(df,'e')
  target <- data.frame(Sample=c('A','B','C'),
                       present=c(0,1,0))
  checkEquals(target,df.result)
  checkEquals(dim(df.result)[1],length(unique(df$Sample)))  
  
} 


test.getAllByPhylum <- function() {
  df <- data.frame(Sample=c('A','A','B','C'),
                   Location = c('u','u','v','v'),
                   phylum=c('tiger','lion','tiger','zebra'),
                   Abundance=c(15,10,20,200),
                   PH_CACL = c(1,1,2,2),
                   PH_H2O = c(1,1,2,2),
                   x = c(10,20,30,40),
                   z = c(11,21,31,41)
                   )
  result <- getAllByPhylum(phylum='zebra',df = df,colnames=c('x','z','PH_H2O','PH_CACL'))
  checkEquals(c('A','B','C'),as.character(result$Sample))
  checkEquals(c(0,0,1),as.numeric(result$present))
  ##
  result <- getAllByPhylum(phylum='tiger',df=df,colnames=c('x','z','PH_H2O','PH_CACL'))
  checkEquals(c('A','B','C'),as.character(result$Sample))
  checkEquals(c(1,1,0),as.numeric(result$present))
  ##
  result <- getAllByPhylum(phylum='lion',df=df,colnames=c('x','z','PH_H2O','PH_CACL'))
  checkEquals(c('A','B','C'),as.character(result$Sample))
  checkEquals(c(1,0,0),as.numeric(result$present))
  checkEquals(c(25,20,200),as.numeric(result$Abundance))
}

test.getGetDFByPhylum() <- function() {
  create.df.species(filepath = filepath,locations = LOCATIONS) %>% 
  getDFByPhylum(phylum='Tardigrada')   %>>% ( ~ temp)
  checkEquals(dim(temp)[1],137)
  checkEquals(dim(temp)[2],17)
  f <- loadData(FILEPATH) 
}

# test.abundance() <- function() {
#   d <- data.frame(Sample=c(1,1),
#                   Species = c('tiger','tiger'),
#                   OTU = c('AA','BB'),
#                   Abundance = c(10,20))
# }

test.getLocationByPhylum <- function() {
  s <- saveLoadedData[]
  ph <- 'Tardigrada'
  #samples <- s[unique(s[which(s$phylum == ph),'Sample']),"Location"]
  
  # s[which(s$Sample %in% samples),"Location"]
  # sum(s[which(s$phylum == ph),'Location'] == 'Lake_Terrasovoe')
  result <- getLocationByPhylum(df=saveLoadedData,phylum='Tardigrada')
  target=c(20,9,3)
  names(target) <- LOCATIONS
  checkEquals(target,result[2,])
}

test.getAbundancesBySample <- function() {
  df <- data.frame(Sample = c('A','A','B','B','A'),
                   Abundance= c(1,2,3,4,5))
  result <- getAbundancesBySample(df)
  target <- data.frame(Sample = c('A','B'),
                       Abundance = c(8,7) )
  checkEquals(target,result)
  
  df <- data.frame(Sample = c('A','A','B','C','B'),
                   Abundance= c(1,2,3,4,5))
  
  result <- getAbundancesBySample(df)
  target <- data.frame(Sample = c('A','B','C'),
                       Abundance = c(3,8,4) )
  checkEquals(target,result)
  
}

test.getLocationByPhylum2  <- function() {
  d <- ps_molten_saved
  d <- d[which(d$Location %in% LOCATIONS),]
  a <- getAbundancesBySample(d)
  low.samples <- a[which(a$Abundance < 1000),'Sample']
  lows <- which(d$Sample %in% low.samples)
  d <- d[-lows,]
  z <- which(d$Abundance == 0)
  d <- d[-z,]

  # t <- which(d$phylum == 'Tardigrada')
  # e <- d[t,]
  # e <- e[!duplicated(e$Sample),]
  # colnames(e)
  # e$present = 1
  # 
  # e2 <- e[,LR.COLUMNS]
  # e2[ e2 == "NA" ] <- NA
  # e3 <- e2[complete.cases(e2),]
  # dim(e3)
  # target <- table(e3$Location)
  # result <- getLocationByPhylum(df=saveLoadedData,phylum='Tardigrada')[2,]
  # checkEquals(as.numeric(target), as.numeric(result))

  for (phyls in c('Chordata','Tardigrada')) {
    this.phyl = phyls
    t <- which(d$phylum == this.phyl)
    e <- d[t,]
    e <- e[!duplicated(e$Sample),]
    colnames(e)
    e$present = 1
    
    e2 <- e[,LR.COLUMNS]
    e2[ e2 == "NA" ] <- NA
    e3 <- e2[complete.cases(e2),]
    dim(e3)
    target <- table(e3$Location)
    result <- getLocationByPhylum(df=saveLoadedData,phylum=this.phyl)[2,]
    #colnames(result) <- LOCATIONS
    checkEquals(as.numeric(target), as.numeric(result))  
  }   
}


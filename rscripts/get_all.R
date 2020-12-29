 
library(glmnet)
library(pROC)

setwd("C:/Users/Michel/Documents/eukaryotes/rscripts")
source('euk_functions.R')
set.seed(123)

#l <- as.data.frame(loadData())
f <- create.df.species()
ph <- sort(unique(f$phylum))

getEst <- function(x,y) {
    m.cv <- cv.glmnet(x=x,y=y,family='binomial',
                      grouped=F,alpha=1,nfolds = length(y))
    best.lambda <- m.cv$lambda.min
    m.lasso <- glmnet(x=x,y=y,lambda = best.lambda, family='binomial',alpha=1)
    return(m.lasso$beta)
}

 



doStats <- function(df,phylum,n.boots=1000,do.boots=F){
  print(phylum)
  beta=NULL
  lasso.ci=NULL
 
  toofew <- FALSE
  
  s <- getAllByPhylum(df=df,phylum=phylum)
  #s <- getAllByPhylumSpecies(df=df,phylum=phylum,species='Mesobiotus furciger')
  
  presence.table <- table(s$distinct.otus>0,s$Location)
  if (sum(presence.table[2,]) < 12 | sum(presence.table[1,]) < 12 )  {
    toofew <- TRUE
    print('too few or too many')
    result <- list(phylum = phylum,
                   presence.table = presence.table,
                   beta=beta,
                   lasso.ci=lasso.ci,
                   xg.acc=xg.acc,
                   auc.test=auc.test,
                   imp=imp)
  } else {
    
    s2 <- s
    ix <- which(colnames(s2) %in% c("Sample", "AMMN","NITR","PHOS","CARB","QUARTZ"))
    s2 <- s2[,-ix]
    s2$Mawson <- 0
    s2[which(s2$Location=="Mawson Escarpment"),"Mawson"] <- 1
    s2$Menzies <- 0
    s2[which(s2$Location=="Mount Menzies"),"Menzies"] <- 1
    s2$Location <- NULL
    s2 <- s2[complete.cases(s2),]
    # presence
    y <- as.numeric(s2$distinct.otus > 0)
    distinct.otus <- s2$distinct.otus
    s2$distinct.otus <- NULL
    # scale all columns
    for ( j in 1:dim(s2)[2]) {
      s2[,j] <- (s2[,j] - mean(s2[,j],na.rm=T)) /sqrt(var(s2[,j],na.rm=T))
    }
    x=as.matrix(s2)
    var.names <- colnames(s2)
    save(var.names,file='varnames.RData')
    n = dim(x)[1]
    p = dim(x)[2]
    sigma = 1
    # choose nfolds=n, to avoid variability from cv folds
    m.cv <- cv.glmnet(x=x,y=y,family='binomial',alpha=1,nfolds = n)
    best.lambda <- m.cv$lambda.min
    best.lambda
    mg <- glmnet(x,y,lambda=best.lambda,family='binomial',alpha=1,thresh = 1e-12)
    beta = coef(mg, x=x, y=y, s=best.lambda, exact=TRUE)
  
    if (do.boots==T) {  
    
    betas.matrix <- matrix(nrow=dim(x)[2],
                       ncol=n.boots)
    
 
    for (b in 1: n.boots) {
      if (b %% 100 == 1 ) print(b)
      ixb <- sample(1:dim(x)[1],replace=T)
      e <- tryCatch( {
        getEst(x[ixb,],y[ixb])
      }, 
      error=function(cond) {
        message(cond)
        return(NA)
      },
      warning=function(cond) {
        #message(cond)
        return(NA)
      })
      if (b==1) {
        bm <- e
      } else {
        bm <- cbind(bm,e)
      }
    }
    # class(bm)
    bm2 <- matrix(bm,nrow=dim(bm)[1],ncol=dim(bm)[2])
    rownames(bm2) <- var.names
    save(bm2,file=paste0('../boots/',phylum,'.RData'))
    
  
   }
    
   
   result <- list(phylum = phylum,
                   presence.table = presence.table,
                   toofew = toofew,
                   beta=beta  
                   )
  }
  return(result)
}


all.stats <- list()
for (i in seq(1 , length(ph))) {
  print(paste(ph[i],Sys.time()))
  all.stats[[i]] <- doStats(df=f,phylum=ph[i],n.boots=1000,do.boots=T)
}
save.image(file='../rdata/allstats.rdata')
load(file='../rdata/allstats.rdata')


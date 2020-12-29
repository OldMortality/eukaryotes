rm(list=ls())


#load(file = 'C:/Users/Michel/Documents/eukaryotes/eukaryotes_nesi/bootsArthropoda.RData')
load(file = 'C:/Users/Michel/Documents/eukaryotes/boots/Arthropoda.RData')
dim(bm2)


library(ggplot2)

colnames <- c("Potassium",	"Sulphur",	"Conductivity"	,"PH",	"RLU",	"Feldspar",	
              "Titanite", "Garnets",
              "Micas","Dolomite",	"Kaolchlor",
              "Calcite","Chlorite",	"Slope",	"Mawson",	"Menzies")



dat <- stack(as.data.frame(t(bm2)))
p1 <- ggplot(dat) + 
  geom_boxplot(aes(y = ind, x = values)) +
  scale_y_discrete(labels=colnames) +
  ggtitle('Arthropods') +
  labs(x = "boostrap estimates of coefficients", y= "")
p1

# proportion of bootstrap samples not equal to zero
z <- apply(bm2, 1, function(c)sum(c!=0))
z <- z/1000  
df.bar <- data.frame(predictor = colnames, prop = z)
p2 <- ggplot(df.bar,aes(y=predictor, x=prop)) + geom_bar(stat="identity",fill='steelblue') +
  labs(x = "Proportion of samples with predictor in the active set", y= "") + 
  ggtitle("Arthropods")
p2


library(matrixStats)
rq <- rowQuantiles(bm2,  probs = c(0.05, 0.95) )
df2 <- data.frame(predictors = colnames,
                  lower = rq[,1],
                  upper = rq[,2])

p2 <- ggplot(df2,aes(y=predictors,yend=predictors,  x=lower, xend=upper)) + 
  geom_segment(stat="identity", col="steelblue", lineend ="round") +
  labs(x = "90% bootstrap confidence intervals") + 
  ggtitle("Arthropods")+
  scale_y_discrete(labels=colnames) + geom_point(col="steelblue") 
p2


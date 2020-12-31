rm(list=ls())

library(ggplot2)
library(stringr) # for str_to_title()
library(gridExtra)

setwd("C:/Users/Michel/Documents/eukaryotes/rscripts")
source('euk_functions.R')
f <- create.df.species()
df.soilData  <-  getSoilData(f,colnames = MYCOLNAMES)


my.color = 'deepskyblue4'

createBoxplot <- function(df,colname,header,ylabel) {
  df[which(df$Location=="Mawson Escarpment"),"Location"] <- "Mawson"
  d <- data.frame(x = as.numeric(df[[colname]]),
                  y = df$Location)
  d <- d[complete.cases(d),]
  p <- ggplot(data=d,aes(y,x)) +
    geom_boxplot() + ggtitle(header) +
    ylab(ylabel) + xlab('')   +
    theme(text = element_text(size=10))
  return(p)
}

headers <-
colnames <- setdiff(MYCOLNAMES,c("Sample", "AMMN","NITR","PHOS","CARB"))


headers <- c('Potassium','Sulphur','Conductivity','PH','RLU',c(str_to_title(colnames[6:length(colnames)])))
ylabels <- c('Colwell  mg/Kg','mg/Kg','dS/M','PH','Relative Light Units',rep('Proportion',9),'?')
# must get unique samples
i = 1
all.plots <- list()
for (i in 1:length(colnames)) {
  p <- createBoxplot(df=f[!duplicated(f$Sample),],
                     colname=colnames[i],
                     header=headers[i],
                     ylabel=ylabels[i])
  p
  all.plots[[i]] <- p
  ggsave(p, file=paste0("../images/boxplots_soildata/",headers[i],'.pdf'), 
         device="pdf",dpi=800,
         width=100,height=70,units='mm')
  
}





# first plot
# page1 <- grid.arrange(all.plots[[1]],
#                            all.plots[[2]],
#                            all.plots[[3]],
#                            all.plots[[4]],
#                            ncol=3)
# ggsave(page1, file=paste0("../images/boxplots_page1"), 
#        device="bmp",dpi=800,
#        width=100,height=70,units='mm')
# 
# # all.plots[[5]],
# # all.plots[[15]],
# 
# 
# # second plot
# page2 <- grid.arrange(all.plots[[6]],
#                       all.plots[[7]],
#                       all.plots[[8]],
#                       all.plots[[9]],
#                       # all.plots[[10]],
#                       # all.plots[[11]],
#                       # all.plots[[12]],
#                       # all.plots[[13]],
#                       # all.plots[[14]],
#                       ncol=2)
# ggsave(page2, file=paste0("../images/boxplots_page2"),
#        device="bmp",dpi=800,
#        width=200,height=100,units='mm')
# 
# # second plot
# page3 <- grid.arrange(all.plots[[10]],
#                       all.plots[[11]],
#                       all.plots[[12]],
#                       all.plots[[13]],
#                       #all.plots[[14]],
#                       ncol=2)
# ggsave(page3, file=paste0("../images/boxplots_page3"),
#        device="bmp",dpi=800,
#        width=200,height=100,units='mm')

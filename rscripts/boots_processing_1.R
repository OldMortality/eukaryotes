
process_bootstraps <- function(phylum,bm2,var.names) {

dim(bm2)
library(ggplot2)
dim(bm2)

not.zero <- vector()
for (i in 1:dim(bm2)[1]) {
  not.zero[i] <- length(which(bm2[i,] != 0))
}
not.zero <- not.zero/ 1000
not.zero

df.prop <- data.frame(x=var.names,
                      y=not.zero)
df.prop[which(df.prop$x=='Menzies'),'x'] <- 'MENZIES'
df.prop[which(df.prop$x=='Mawson'),'x'] <- 'MAWSON'
head(df.prop)
# get the plot ordered alphabetically from top to bottom
df.prop$x <- factor(df.prop$x)
df.prop$x <- factor(df.prop$x,ordered=T,levels=rev(levels(df.prop$x)))

my.color = 'deepskyblue4'
#reorder(position, desc(position))
plot.zeros <- ggplot(data=df.prop) + 
  geom_col(aes(x=x,y=y),fill=my.color, colour=my.color) +
  coord_flip()   +
  geom_hline(yintercept=0.95,colour='red',linetype = "dashed") +
  ylab('Percentage of samples not zero') +
  xlab('') + 
  expand_limits(y=c(0,1)) + 
  labs(title = phylum, subtitle = "bootstrap probability of not 0")
plot.zeros

lower <- vector()
upper <- vector()
for (i in 1:dim(bm2)[1]) {
  ci <- quantile(bm2[i,],c(0.025,0.975),na.rm=T)
  lower[i] <- ci[1]
  upper[i] <- ci[2]
}
d <- data.frame(varnames = var.names,lower=lower,upper=upper)

#get the plot ordered alphabetically from top to bottom
d$varnames[which(d$varnames=='Menzies')] <- 'MENZIES'
d$varnames[which(d$varnames=='Mawson')] <- 'MAWSON'
d$varnames <- factor(d$varnames)
d$varnames <- factor(d$varnames,ordered=T,levels=rev(levels(d$varnames)))

plot.ci <- ggplot(data=d) + geom_segment(aes(y=varnames,
                              yend=varnames,
                              x=lower,
                              xend=upper), colour=my.color) +
  ylab('') +
  xlab('95% confidence interval') +
  geom_vline(xintercept=0,colour='red',linetype = "dashed") +
  
  labs(title = "", subtitle = "") 
plot.ci

cis <- data.frame(x=var.names,
                 lower=round(lower,3),
                 upper=round(upper,3))  




return(list(phylum = phylum,
            plot.zeros = plot.zeros,
            plot.ci = plot.ci,
            cis = cis,
            df.prop = df.prop))

}


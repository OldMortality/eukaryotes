
process_bootstraps <- function(phylum,bm2,var.names) {

dim(bm2)
library(ggplot2)
dim(bm2)

# dat <- stack(as.data.frame(t(bm2)))
# colnames(dat)[1] <- 'Estimates'
# plot.est <-ggplot(dat) + 
#   geom_boxplot(aes(x = ind, y = Estimates)) +
#   coord_flip() + 
#   scale_x_discrete(name='predictors',
#                    labels=var.names) +
#   ylab('LASSO estimate') +
#   xlab('predictor') +
#   labs(title = phylum, subtitle = "Bootstrap estimates")


not.zero <- vector()
for (i in 1:dim(bm2)[1]) {
  not.zero[i] <- length(which(bm2[i,] != 0))
}
not.zero <- not.zero/ 1000
not.zero

df.prop <- data.frame(x=var.names,
                      y=not.zero)
head(df.prop)

my.color = 'deepskyblue4'

plot.zeros <- ggplot(data=df.prop) + geom_col(aes(x=x,y=y),fill=my.color, colour=my.color) +
  coord_flip() + 
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
            cis = cis))

}


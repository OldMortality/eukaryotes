library(binaryLogic)
b <- as.character(as.binary(8))
which(b=="1")


new_ph <- function(x,y,p) {
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


N <- 3 # number of phyla
p <- c(0.2,0.1,0.7)
stopifnot(length(p)==N)
T <- matrix(0,nrow=2^N,ncol=2^N)
for (i in seq(1,2^N-1)) {
  T[i,i] <- sum(p * as.numeric(as.character(as.binary(i-1,n = N)))) 
}
T[2^N,2^N] <- 1

for (i in seq(1,2^N-1)) {
  for (j in seq(i + 1,2^N)) {
    if (i != j) {
      T[j,i] <- new_ph(x=i,y=j,p=p)  
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
t_vec <- M %*% one_vector
# expected number of steps
t_vec
# variance of the expected number of steps
(2 *  M - I ) %*% t_vec - t_vec * t_vec




  ns <- vector()  
  for (sim in 1:10000) {
    seen <- matrix(0,nrow=N,ncol=1)
    for (i in 1:100) {
      y <- rmultinom(1,size=1,prob=p)
      seen <- y + seen
      for(j in 1:N) {
        seen[j] <- min(seen[j],1)
      }
      if (identical(seen,matrix(1,nrow=N,ncol=1))) {
        ns[sim] <- i
        break
      }
    }
  }
quantile(ns,0.95)  
mean(ns)

  
d <- data.frame(x=as.numeric(ns))
hist(d$x,probability = T)
var(ns)
sd(ns)
mean(ns)
mean(ns) + 2 * sd(ns)


library(ggplot2)

c <- data.frame(x=rchisq(1000,df=mean(d$x)))

#
ggplot() + geom_density(data=c,aes(x=x),col='red') +
geom_density(data=d,aes(x=x))

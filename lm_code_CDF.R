library(distrEx)
library(mvtnorm)

set.seed(111) ## not the original seed used in paper
#### data ####
bet0 <- c(1,0.5,0.25,-1,0.75,-1.5,2,1.6,1,0.8)

p <- 10     ## no. of predictors
n <- 5000  ## training sample size
m <- 50    ## test sample size
x <- matrix(runif(n*p), n, p)  
x.test <- matrix(runif(m*p), m, p)

fsimu <- function(x,beta) return(x%*%beta)
y <- apply(x, 1, fsimu, beta=bet0) + rnorm(n,0,1) 
obs <- scale(y)
sy <- attr(obs, "scaled:scale")
my <- attr(obs, "scaled:center")

nstr  <- 500      ## number of batches
Xlist <- list()
ylist <- list()


ns <- n/nstr      ## obs per batch
for(i in 1:nstr){
  Xlist[[i]]<-x[((i-1)*ns+1):(i*ns),]  ## predictors in the ith batch
  ylist[[i]]<-y[((i-1)*ns+1):(i*ns)]   ## responses in the ith batch
}

#posterior samples 1,2,...t
nsamp <- 500
str_sig2.store <- replicate(nstr,list(rep(NA,nsamp)))
str_beta.store <- replicate(nstr,list(matrix(NA,nsamp,p)))
cdfstr_sig2.store <- replicate(nstr,list(rep(NA,nsamp)))
cdfstr_beta.store <- replicate(nstr,list(matrix(NA,nsamp,p)))

### c-df
yy  <- 0              #streaming: y'y
XX  <- matrix(0,p,p)  #streaming: X'X
XY  <- matrix(0,p,1)  #streaming: X'y
S1xx <- matrix(0,p,p)
S1xy <- matrix(0,p,1)
S2xx <- 0
S2xy <- 0

#prior initializations
Sigmabet <- diag(p)
ak <- 0
bk <- 0
b <- 2

for(j in 1:nstr){
  x1 <- Xlist[[j]]
  y1 <- ylist[[j]]

  ## update quantities
  #yy  <- yy + crossprod(y1,y1)  #[1], cumul store y'y
  #XX  <- XX + crossprod(x1,x1)  #[p x p], cumul store  X'X
  #XY  <- XY + crossprod(x1,y1)  #[p x 1], cumul store  X'y
  
  ## initialize at the last iterate of the previous time point
  if(j==1) {
        sig22 <- 1
        beta <-  bet0
  } else{
        sig22 <- cdfstr_sig2.store[[j-1]][nsamp] 
        beta <- cdfstr_beta.store[[j-1]][nsamp,] 
  }


  ## params
  ak <- ak + length(y1)/2

  S1xx <- S1xx + crossprod(x1,x1)/sig22   ## Surrogate CSS
  S1xy <- S1xy + crossprod(x1,y1)/sig22   ## Surrogate CSS

  vari1 <- chol2inv(chol(S1xx+diag(p)))
  mu1 <- vari1 %*% S1xy


  ## sample
  for(i in 1:nsamp){
    cdfstr_beta.store[[j]][i,] <- c(rmvnorm(1,mu1,vari1))  ## beta update
  }

  mymu <- colMeans(cdfstr_beta.store[[j]]) ## beta point estimate
  S2xx <- c(crossprod(y1-x1%*%mymu))+S2xx  ## Surrogate CSS
  bk1 <- b+S2xx/2   

  for(i in 1:nsamp){
    cdfstr_sig2.store[[j]][i] <- 1/rgamma(1,ak,bk1)  ## sigma update
  }

  cat('done stream #',j,"...\n")
 }

### Sequential Markov Chain Monte Carlo
yy  <- 0              # Sufficient Statistics: y'y
XX  <- matrix(0,p,p)  # Sufficient Statistics: X'X
XY  <- matrix(0,p,1)  # Sufficient Statistics: X'y
S1xx <- matrix(0,p,p)
S1xy <- matrix(0,p,1)
S2xx <- 0
S2xy <- 0

#prior initializations
Sigmabet <- diag(p)
ak <- 0
bk <- 0
b <- 2

for(j in 1:nstr){
  x1 <- Xlist[[j]]
  y1 <- ylist[[j]]

  ## update quantities
  yy  <- yy + crossprod(y1,y1)  #[1], cumul SS store
  XX  <- XX + crossprod(x1,x1)  #[p x p], cumul SS store
  XY  <- XY + crossprod(x1,y1)  #[p x 1], cumul SS store

  ## initialize at the last iterate of the previous time point
  if(j==1) {
        sig22 <- 1
        beta <-  bet0
  } else{
        sig22 <- cdfstr_sig2.store[[j-1]][nsamp]
        beta <- cdfstr_beta.store[[j-1]][nsamp,]
  }


  ## params
  ak <- ak + length(y1)/2

  ## sample
  for(i in 1:nsamp){
    bk <- b+(yy-2*crossprod(beta,XY)+t(beta)%*%XX%*%beta)/2
    sig2 <- 1/rgamma(1,ak,bk)      ## sigma update
    str_sig2.store[[j]][i] <- sig2 ## sigma stored

    vari <- chol2inv(chol(XX/sig2+diag(p)))
    mu <- vari%*%(XY/sig2)
    beta <- c(rmvnorm(1,mu,vari))   ## beta update
    str_beta.store[[j]][i,] <- beta ## sigma stored
  }
  
 cat('done stream #',j,"...\n")
}


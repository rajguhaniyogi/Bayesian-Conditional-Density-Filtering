# library(mvtnorm)
library(MASS)
# library(BayesBridge)
library(msm)
library(truncnorm)

### defs
ff <- function(xx, yy, bet) {
    val <- rep(0, length(yy))    
    means  <- xx %*% bet
    numer  <- (2 * yy - 1) * dnorm(-means)
    denom  <- sapply(1:nn, function(ii) pnorm(q = -means[ii], lower.tail = 1 - yy[ii]))    
    val <- means + numer / denom
    
    return(val)
}

### define dataset 
simulate <- function(N, bet) {
    p <- length(bet)
    X <- matrix(rnorm(N * p, sd = 0.2), N, p)
    ystar <- X %*% bet + rnorm(N)
    y <- as.numeric(ystar > 0)
    print(table(y))
    
    list(y=y, X=X, ystar=ystar)
}

#N <- length(y.new)
#p <- ncol(X.new)

N <- 10e3
p <- 500
beta.tr <- c(7,-7,-4,4,-3,3,-3,3,-2,2) / 2
beta.tr <- c(beta.tr, runif(p - 10, min = -0.75, max = 0.75))

D <- simulate(N, beta.tr)  ## Simulated data

### init & run
nstr <- 100    ## number of batches
n <- N / nstr  ## sample size per batch
#budget <- max(0.10 * N, round(p * log(p) / 100) * 100)
budget <- 0.2*N ## size of the moving window
myitr <- rep(500, nstr) ## number of gibbs draws for c-df

lo <- c(-Inf, 0)
hi <- c(0, Inf)

out.id <- 1
css <- rep(0, p) 
ss  <- matrix(0, p, p) 

css.store <- array(0, dim = c(nstr - (budget / n), p))
beta.store <- array(0, dim = c(p, nstr, max(myitr)))
bet  <- rep(0, p)
Z    <- rnorm(n, 0, 1)

xx <- yy <- list()
for(i in 1:nstr){
    xx[[i]] <- D$X[((i-1)*n+1):(i*n),]  ## predictors in the i-th batch
    yy[[i]] <- D$y[((i-1)*n+1):(i*n)]   ## response in the i-th batch
}

### C-DF
for(i in 1:nstr){

    updateIds  <- 1:i
    budget.str <- floor(budget / n)
    if(i > budget.str) {
        updateIds <- seq(i-budget.str+1, i, by=1)
        
        out.id <- (i-budget.str)        
        css <- css + crossprod(xx[[out.id]], Z[1:n]) ## update surrogate CSS
        css.store[out.id, ] <- crossprod(xx[[out.id]], Z[1:n])
    }
    
    ii <- seq((updateIds[1]-1) * n + 1, i * n, by = 1)
    myX <- D$X[ii, ]
    myy <- D$y[ii]
    nn  <- length(myy)
    
    ## update stats
    ss <- ss + crossprod(xx[[i]], xx[[i]])  ## surrogate CSS

    ## sample z
    post.var <- chol2inv(chol(ss + diag(p)))
    chol.var <- chol(post.var)
    
    LO  <- lo[myy + 1]
    HI  <- hi[myy + 1]

    myZ <- rep(0, nn)
    mybet <- rep(0, p)
    itr <- myitr[i]
    
    if(out.id > 1) {
        aa <- 1
        css <- apply(css.store[1:out.id, ] * aa, 2, sum)
    } else {
        aa <- 1
        css <- css.store[1:out.id, ]
    }    

    ## gibbs 1...itr
    for(it in 1:itr) {
        ## sample z
        Z <- rtruncnorm(1, LO, HI, myX %*% bet, 1)       ## update latent variable 
        if(sum(is.nan(Z))) cat('NAs produced\n')
        myZ <- myZ + (Z * ifelse(it > 0.5 * itr, 1, 0))   
        
        ## sample beta
        mycss <- css + crossprod(myX, Z)
        bet <- as.numeric(crossprod(post.var, mycss)) + chol.var %*% rnorm(p)  ## beta update
        beta.store[, i, it] <- bet
    }
    
    Z <- myZ / (0.5 * itr)
    cat(i, "\t", "ngibbs: ", it, "\n")
    cat(round(apply(beta.store[1:10, i, ], 1, mean), 2), "\n")
}


### Full Gibbs sampler
ss <- matrix(0, p, p)
strs <- c(1:nstr)
bet.gibbs.store <- array(0, dim = c(p, length(strs), max(myitr)))
bet  <- rep(0, p)
Z    <- rnorm(n, 0, 1)

for(i in seq_along(strs)) {
    itr <- myitr[i]
    ii <- seq(1, strs[i] * n, by = 1)
    myX <- D$X[ii, ]
    myy <- D$y[ii]
    nn  <- length(myy)
    
    ss <- crossprod(myX, myX)
    post.var <- chol2inv(chol(ss + diag(p)))
    chol.var <- chol(post.var)
    LO  <- lo[myy + 1]
    HI  <- hi[myy + 1]
    
    for(it in 1:itr) {
        Z <- rtnorm(n=nn, mean=c(myX %*% bet), sd=1, lower=LO, upper=HI)  ## update latent variables
        bet <- crossprod(post.var, t(myX) %*% Z) + chol.var %*% rnorm(p)  ## update beta
        
        ## store
        bet.gibbs.store[, i, it] <- bet
        
    }
    
    cat(i, "\t", "ngibbs: ", it, "\n")
    cat(round(apply(bet.gibbs.store[1:10, i, ], 1, mean), 2), "\n")
}


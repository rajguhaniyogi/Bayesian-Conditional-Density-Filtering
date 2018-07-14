library(distrEx)
library(pscl)
set.seed(2681)

### simulated data
nmax<-1000         ## number of data batches

it <- 500          ## number of iterations per batch
true.mu <- 4       
true.sigma <- 100  
true.tau <- 0.01
n <- 100           ## data points per batch
k <- 10            ## number of local effects 

a.m <- 2
b.m <- 1

zeta.true <- rnorm(k,true.mu,sqrt(true.tau))  ## true local effects
y <- rnorm(n*nmax*k,kronecker(zeta.true,rep(1,n*nmax)),sqrt(true.sigma))  ## generated response
mat.y <- matrix(y,n*nmax,k) 


### Conditional Density Filtering
css3 <- rep(0,k)
css4 <- 0

testsse.cdf <- rep(0, nmax)
mu.cdf <- matrix(NA,it,nmax)
tau.cdf <- matrix(NA,it,nmax)
zeta.cdf <- list()
sigma.cdf <- matrix(NA,it,nmax)

for(h in 1:nmax){
	suby <- mat.y[(n*(h-1)+1):(n*h),]  ## data in the h-th batch
	zeta.cdf[[h]] <- matrix(NA,it,k)   
      zet1 <- zeta.cdf[[h]]
    
    if(h==1) {
        mymu <- 0
        mytau <- 0.01
        myzet <- rep(0, k)
        mysig <- 10
        
        tau.cdf[1,h] <- mytau
        zet1[1,] <- myzet
        mu.cdf[1,h] <- mymu
        sigma.cdf[1,h] <- mysig       
    } else{
        tau.cdf[1,h] <- mytau
        zet1[1,] <- zeta.cdf[[h-1]][it, ]
        mu.cdf[1,h] <- mymu
        sigma.cdf[1,h] <- sigma.cdf[it, h-1]
    }
      ## Update zeta and sigma^2 under C-DF 
	for(i in 2:it){
		for(j in 1:k){
		mean.zetaj1 <- (css3[j] + sigma.cdf[i-1,h] * mymu) / (n*h*mytau + sigma.cdf[i-1,h])
		var.zetaj1 <- (mytau * sigma.cdf[i-1,h]) / (n*h*mytau + sigma.cdf[i-1,h])
		zet1[i,j] <- rnorm(1, mean.zetaj1, sqrt(var.zetaj1)) ## zeta update
		
	}
	
	mycss2 <- css4 + sum(colSums(suby^2)) - 2*sum(colSums(suby) * zet1[i, ]) + n*sum(zet1[i, ]^2)
	sigma.cdf[i,h] <- 1 / rgamma(1, n*h*k / 2, mycss2 / 2)  ## sigma^2 update

}

    myzet <- apply(zet1, 2, mean)
    mysig <- mean(sigma.cdf[, h]) ##not used in css
    
    ## Update level 2 params mu and tau^2
    for(i in 2:it) {
        
        mu.cdf[i,h] <- rnorm(1, mean(myzet), sqrt(tau.cdf[i-1,h] / k))
        tau.cdf[i,h] <- 1 / rgamma(1, a.m + k/2, b.m + sum((myzet - mu.cdf[i,h])^2) / 2)
        
        a_tau <- a.m+k/2
        b_tau <- b.m + sum((myzet - mu.cdf[i,h])^2)
    }
    
    ## finalize some params 
    mymu  <- mean(mu.cdf[, h]) ## point estimate of mu 
    mytau <- mean(tau.cdf[, h]) ## point estimate of tau^2
    
    css3 <- css3 + mytau * colSums(suby)  ## Surrogate CSS
    css4 <- css4 + sum(colSums(suby^2)) -2*sum(colSums(suby) * myzet) + n*sum(myzet^2)  ## Surrogate CSS
    
    ## store dist by batch
    zeta.cdf[[h]] <- zet1

    cat('sim:::::::::::::::::::: ',h,"\n")
}

### Sequential Markov Chain Monte Carlo
css <- rep(0,k)
css2 <- 0

mu <- matrix(NA,it,nmax)
tau <- matrix(NA,it,nmax)
zeta <- list()
sigma <- matrix(NA,it,nmax)

for(h in 1:nmax){
	suby <- mat.y[(n*(h-1)+1):(n*h),] ## data at h-th batch
	css <- css+colSums(suby) ## Sufficient Statistic
	css2 <- css2+sum(colSums(suby^2)) ## Sufficient Statistic
	
	## MCMC 
	zeta[[h]] <- matrix(NA,it,k) 
	zet <- zeta[[h]]
	mu[1,h] <- 1
	zet[1,] <- rnorm(k,0,1)
	sigma[1,h] <- 100
	
	if(h==1) {
		tau[1,h] <- .01
		zet[1,] <- rnorm(k,0,1)
		mu[1,h] <- 0
		sigma[1,h] <- 100
	} else{
		tau[1,h] <- tau[it,h-1]
		old.zet <- zeta[[h-1]]
		zet[1,] <- old.zet[it,]
		mu[1,h] <- mu[it,h-1]
		sigma[1,h] <- sigma[it,h-1]
	}
	
    ## Update parameters in SMCMC
    for(i in 2:it) {
        for(j in 1:k){

        mean.zetaj <- (tau[i-1,h]*css[j]+sigma[i-1,h]*mu[i-1,h])/(n*h*tau[i-1,h]+sigma[i-1,h])
        var.zetaj <- (tau[i-1,h]*sigma[i-1,h])/(n*h*tau[i-1,h]+sigma[i-1,h])
        zet[i,j] <- rnorm(1,mean.zetaj,sqrt(var.zetaj))  ## zeta update
        }
       
        mu[i,h] <- rnorm(1,mean(zet[i,]),sqrt(tau[i-1,h]/k))           ## mu update
        tau[i,h] <- rigamma(1,a.m+k/2,b.m+sum((zet[i,]-mu[i,h])^2)/2)  ## tau^2 update
        b_tau1 <- b.m+sum((zet[i,]-mu[i,h])^2)
                
        a_sig <- n*h*k
        b_sig <- css2 -2*sum(css * zet[i,]) + n*h*sum(zet[i,]^2)
        sigma[i,h] <- 1/rgamma(1,a_sig/2,b_sig/2)                      ## sigma^2 update
    }
    
    ## finalize some params 
    	zeta[[h]] <- zet

	cat('sim:::::::::::::::::::: ',h,"\n")
}

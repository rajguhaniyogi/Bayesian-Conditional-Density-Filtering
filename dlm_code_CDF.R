rm(list = ls())
library(pscl)
library(mvtnorm)
library(msm)


## functions to help compute "scss" for C_DF
ff <- function(phi, thet) { 
	##only when length(thet) > 2
	t1 <- thet[-length(thet)]
	t2 <- thet[-1]

	return(sum((t2 - phi * t1)^2))	
}

ff2 <- function(thet) { 
	##only when length(thet) > 2
	t1 <- thet[-length(thet)]
	t2 <- thet[-1]

	return(sum(t2 * t1))
}

### data generation
n.pred   <- 500   
n.data   <- 5000
n.split  <- 50
n.per    <- 100 # n.data / n.split

tau.true <- 2
phi.true <- 0.7
sigma.true <- .01
theta.init <- 3 ##theta0
set.seed(999)

theta.true <- numeric()
y <- numeric()
y.pred <- list()

for(i in 1:n.data) {
	theta.int <- theta.init
	theta.true[i] <- rnorm(1,phi.true*theta.int,sqrt(tau.true))
	
	if( i%%n.per == 0) {
  		z <- rnorm(n.pred+1,theta.true[i],sqrt(sigma.true))
    	y[i] <- z[1]
    	y.pred[[(i/n.per)]] <- z[-1]
 	} else {
		y[i] <- rnorm(1,theta.true[i],sqrt(sigma.true))
	}
  
  	theta.init <- theta.true[i]
}

### C-DF
a0 <- 2
b0 <- 2
N.particles <- 1000
n.gibbs <- N.particles/2   ## Gibbs samples
N.MH    <- N.particles * 5 ## Metropolis Hastings samples

theta   <- rep(0, n.per)
tau2 <- 1
sig2 <- 0.1
phi  <- 0.6

scss1 <- 0 ## sum of: (y_t - theta_t)^2 
scss2 <- 0 ## sum of: (theta_{t+1} - phi * theta_t)^2
scss3 <- 0 ## sum of: theta_t^2
scss4 <- 0 ## sum of: theta_t * theta_{t+1}

tau2.vec <- numeric()
sig2.vec <- numeric()
phi.vec  <- numeric()

theta.store <- NULL
tau2.store  <- NULL
sig2.store  <- NULL
phi.store   <- NULL

## loop
idx <- c(1)
beg <- Sys.time()
theta.hat <- c()

beg <- Sys.time()

for(i in 1:n.data) {

	## get idx 
	if(i <= n.per) {
		idx <- c(1:i)
	} else {
		idx <- seq(from = i - n.per + 1, to = i, by = 1)
		theta.out <- mean(theta.mat[, 1]) ## outgoing theta fix
    	      theta.hat <- c(theta.hat, theta.out)
		
		## update surrogate CSS
		scss1 <- scss1 + (y[i-n.per] - theta.out)^2
		scss2 <- ifelse(length(theta.hat) >= 2, ff(phi, theta.hat), 0)
		scss3 <- scss3 + theta.out^2	
		scss4 <- ifelse(length(theta.hat) >= 2, ff2(theta.hat), 0)
		
		## store parameters
		if(i %% 100 == 0) {
			theta.store <- rbind(theta.store, theta.mat[, 1]) ## store theta
			tau2.store  <- rbind(tau2.store, tau2.vec) ## store tau2
			sig2.store  <- rbind(sig2.store, sig2.vec) ## store sig2
			phi.store   <- rbind(phi.store, phi.vec)   ## store phi
		}
	}
	
	## do C-DF sampling
	theta.mat <- NULL
	
	for(s in 1:n.gibbs) {
		
		## sample theta		
		for(j in seq(length(idx), to = 1, by = -1)) {
			myid <- idx[j]
			
            ## care with index of theta
			if(j == length(idx)) {
                    if(i == 1) {
                      theta.prev <- theta.init
                    } else {
                    theta.prev <- ifelse(s == 1, theta[j], theta[j-1])
                    }
			    
			  omega  <- 1 / (1 / sig2 + 1 / tau2)
			  numer  <- y[myid] / sig2 + phi * theta.prev / tau2
			} else {
			    if(j == 1) {
			        theta.prev <- ifelse(i <= n.per, theta.init, theta.out)
			    } else {
			        theta.prev <- ifelse(s == 1, theta[j], theta[j-1])
			    }
			    
			    omega  <- 1 / (1 / sig2 + (phi^2 + 1) / tau2)
			    numer  <- y[myid] / sig2 + phi * (theta.prev + theta[j + 1]) / tau2
			  }
			
			theta[j] <- rnorm(1, numer * omega, sqrt(omega))
		}
        
            qq <- ifelse(i > 1, ff(phi, theta), 0)
		
		## sample sig2
		a.sig2 <- a0 + i / 2
		b.sig2 <- b0 + (scss1 + sum((y[idx] - theta[theta!=0])^2)) / 2
		sig2    <- 1 / rgamma(1, a.sig2, b.sig2)
            		
		## sample tau2
		a.tau2 <- a0 + i / 2
		b.tau2 <- b0 + (scss2 + qq) / 2 ## see qq calculated above
		tau2    <- 1 / rgamma(1, a.tau2, b.tau2)
		
		## sample phi
		qq2 <- ifelse(i > n.per, theta.out * theta[1], 0)
		my.sum2.theta <- scss3 + sum(theta^2)
		my.sumX.theta <- scss4 + qq2 + ff2(theta)
		
		prop.phi <- rtnorm(N.MH, phi, 1, -1, 1)
		log.accept.prob.num <- dnorm(prop.phi, my.sumX.theta / my.sum2.theta, sqrt(tau2 / my.sum2.theta), log=T)
		log.accept.prob.den <- dnorm(phi, my.sumX.theta / my.sum2.theta, sqrt(tau2 / my.sum2.theta), log=T)
		
		prob.to.accept <- log(runif(N.MH))
		accept.prob    <- log.accept.prob.num - log.accept.prob.den
		accepted.phi   <- prop.phi[which(prob.to.accept < accept.prob)]
		
		if(length(accepted.phi) > 0){
		    phi <- mean(accepted.phi)
		}
		
		## store updates
		theta.mat   <- rbind(theta.mat, theta)
		sig2.vec[s] <- sig2
		tau2.vec[s] <- tau2
		phi.vec[s]  <- phi
	}
	
	## print
	if(i %% 5 == 0) {
		cat('it: ', i, '\t tau2: ', tau2, '\t sigma2: ', sig2, '\t phi: ', phi, '\n')
	}
	
}

time1 <- Sys.time()-beg
      

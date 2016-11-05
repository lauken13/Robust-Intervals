library(LaplacesDemon,quietly=TRUE,warn.conflicts = FALSE,verbose=FALSE)
library(rjags,quietly=TRUE,warn.conflicts = FALSE,verbose=FALSE)
library(ggplot2,quietly=TRUE,warn.conflicts = FALSE,verbose=FALSE)
library(sn,quietly=TRUE,warn.conflicts = FALSE,verbose=FALSE)
library(weights,quietly=TRUE,warn.conflicts = FALSE,verbose=FALSE)
library(dplyr)
library(HDInterval)
###################################################################################
#                           FUNCTIONS
###################################################################################


createData<-function(sample.size,gval,hval,contProp,contSpread=NA,error){
  if(error=='normal'){
    x<-c(ghdist(sample.size*(1-contProp),gval,hval),rnorm(sample.size*contProp,0,contSpread))
  }else if(error=='normal.bias'){
    x<-c(ghdist(sample.size*(1-contProp),gval,hval),rnorm(sample.size*contProp,10,1))
  }
    return(x)
}

#### model 1 ####

ghdist<-function (n, g = 0, h = 0) 
{
  x <- rnorm(n)
  if (g > 0) {
    ghdist <- (exp(g * x) - 1) * exp(h * x^2/2)/g
  }
  if (g == 0) 
    ghdist <- x * exp(h * x^2/2)
  ghdist
}

BBmean <- function(x,nsamp=1000) {
  samples <- BayesianBootstrap(x,n=nsamp,Method=weighted.mean)
  samples <- unlist(samples)
  return(samples)
}


#### model 2 ####

BBtrim <- function(x,nsamp=1000) {
  
  weighted.trimmed.mean <- function(x,w,trimmed=.2){
    ord <- order(x)
    x <- x[ord]
    w <- w[ord]
    w <- truncate.weights(w,trimmed)
    m <- weighted.mean(x,w)
    return(m)
  }
  
  truncate.weights <- function(w,trimmed=.2){
    tw <- w
    lsum <- 0
    ind <- 0
    while(lsum < trimmed) {
      ind <- ind+1
      lsum <- lsum + w[ind] 
      tw[ind] <- max(0,lsum-trimmed)
    }
    rsum <- 0
    ind <- length(w)+1
    while(rsum < trimmed) {
      ind <- ind-1
      rsum <- rsum + w[ind] 
      tw[ind] <- max(0,rsum-trimmed)
    }  
    tw <- tw/sum(tw)
    return(tw)
  }
  
  samples <- BayesianBootstrap(x,n=nsamp,Method=weighted.trimmed.mean)
  samples <- unlist(samples)
  return(samples)
}

#### model 3 ####

normalModel <- function(x,nsamp=1000,thin=5,burnin=1000) {
  
  jagsModelString <- 
    "model {
  tau ~ dgamma(.0001,.0001) # prior over the precision
  mu ~ dnorm(0,.0001) # prior over mean
  for( i in 1:N) {
  X[i] ~ dnorm(mu,tau)
  }
}"
  
  jagsData <- list(
    X = x-mean(x),
    N = length(x)
  )
  
  jagsModel <- jags.model(
    file = textConnection(jagsModelString), 
    data = jagsData, 
    n.adapt = burnin,
    quiet=TRUE 
  )
  
  samples <- jags.samples(
    model = jagsModel, 
    variable.names = "mu", 
    n.iter = nsamp*thin,
    thin = thin,
    progress.bar="none"
  )
  
  mu <- as.matrix( samples[[1]] )
  samples <- mu[,1] + mean(x)
  return(samples)
  }

#### model 4 ####

CN.mean <- function(x,nsamp=1000,thin=5,burnin=1000) {
  jagsModelString <- 
    "model {
  tau ~ dgamma(.0001,.0001) # prior over the precision
  sig <- 1/sqrt(tau) # standard deviation
  mu ~ dnorm(0,.0001) # prior over mean
  htmp ~ dexp(.1) # the scale: contaminants expected to have 10x width
  h <- htmp + .01 # SIGH 
  p ~ dbeta(1,9) # contaminant probability: 10% contaminant
  for( i in 1:N) {
  z[i] ~ dbern(p) # is contaminant?
  sigma[i] <- (1-z[i])*sig + z[i]*h*sig
  precision[i] <- 1/sigma[i]^2
  X[i] ~ dnorm(mu,precision[i]) 
  }
}"
  
  jagsData <- list(
    X = x-mean(x),
    N = length(x)
  )
  
  jagsModel <- jags.model(
    file = textConnection(jagsModelString), 
    data = jagsData, 
    n.adapt = burnin,
    quiet=TRUE 
  )
  
  samples <- jags.samples(
    model = jagsModel, 
    variable.names = c("mu","z"), 
    n.iter = nsamp*thin,
    thin = thin,
    progress.bar="none"
  )
  
  mu <- as.matrix( samples[[1]])
  samples.mu <- mu[,1] + mean(x)
  #samples.cont<-samples[[2]]
  
  estimates<-apply(samples[[2]],1,mean)
  no.cont<-length(estimates[estimates>.5])
  
  #pop.mean<-mean(samples.mu)
  #ci.low <- quantile(samples.mu,probs=.025)
  #ci.high <- quantile(samples.mu,probs=.975)
  Out<- list(posterior=samples.mu, number.contaminants=no.cont)
  return(Out)
  }

#### model 5 ####

tModel <- function(x,nsamp=1000,thin=5,burnin=1000) {
  
  jagsModelString <- 
    "model {
  tau ~ dgamma(.0001,.0001) # prior over the precision
  mu ~ dnorm(0,.0001) # prior over mean
  df ~ dexp(1) #prior over the degrees of freedom
  k<-df+.001
  for( i in 1:N) {
  X[i]~dt(mu,tau,k) 
  }
}"
  
  jagsData <- list(
    X = x-mean(x),
    N = length(x)
  )
  
  jagsModel <- jags.model(
    file = textConnection(jagsModelString), 
    data = jagsData, 
    n.adapt = burnin,
    quiet=TRUE 
  )
  
  samples <- jags.samples(
    model = jagsModel, 
    variable.names = "mu", 
    n.iter = nsamp*thin,
    thin = thin,
    progress.bar="none"
  )
  
  mu <- as.matrix( samples[[1]] )
  samples <- mu[,1] + mean(x)
  return(samples)
  }


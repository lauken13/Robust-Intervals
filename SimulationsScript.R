library(LaplacesDemon,quietly=TRUE,warn.conflicts = FALSE,verbose=FALSE)
library(rjags,quietly=TRUE,warn.conflicts = FALSE,verbose=FALSE)
library(ggplot2,quietly=TRUE,warn.conflicts = FALSE,verbose=FALSE)
library(sn,quietly=TRUE,warn.conflicts = FALSE,verbose=FALSE)
library(weights,quietly=TRUE,warn.conflicts = FALSE,verbose=FALSE)
library(dplyr)

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
  samples.cont<-samples[[2]]
  
  estimates<-apply(samples[[2]],1,mean)
  no.cont<-length(estimates[estimates>.5])
  
  pop.mean<-mean(samples.mu)
  ci.low <- quantile(samples.mu,probs=.025)
  ci.high <- quantile(samples.mu,probs=.975)
  Out<- list(population.mean=pop.mean, credible.interval=c(ci.low,ci.high),number.contaminants=no.cont)
  return(Out)
  }

###################ERROR MODEL:##############################################
###################################################################################
                                #ACTUAL SCRIPT
####################################################################################
#df<-expand.grid(Iteration=1:100, correct='AAA',width=-1,Gval=c(0,0.5),Hval=c(0,0.5),contProp=c(0,.1,.2,.3,.4,.5),contSpread=c(2,5,10),sample.size=c(10,25,50,100,250),error=c('normal','normal.bias','unif','unif.bias','dp'),model=c('percentile_t','t_test'))
#df<-expand.grid(Iteration=1:100, correct='AAA',width=-1,Gval=0,Hval=0,contProp=c(0,.1,.2,.3,.4,.5),contSpread=2,sample.size=c(25,50,100,250),error=c('normal.bias','normal','dp'),model=c('BB mean','BB trim','normal','cont','percentile_t','t_test'))
df<-expand.grid(Iteration=1:1000,seeds=-1,contSkew='AAA',BBmeanCorrect='AAA',BBtrimCorrect='AAA',normalCorrect='AAA',contCorrect='AAA',BBmeanWidth=-1,BBtrimWidth=-1,normalWidth=-1,contWidth=-1,no.contaminants=-1,Gval=0,Hval=0,contProp=c(0,.1,.2,.3,.4,.5),contSpread=10,sample.size=c(20,50,100,250),error=c('normal.bias','normal'))
rand.seed<-sample.int(size=nrow(df),n=100000,replace=FALSE)
df$BBmeanCorrect<-as.character(df$BBmeanCorrect)
df$BBtrimCorrect<-as.character(df$BBtrimCorrect)
df$normalCorrect<-as.character(df$normalCorrect)
df$contSkew<-as.character(df$contSkew)
df$contCorrect<-as.character(df$contCorrect)
df$seeds<-rand.seed
for (r in 1:nrow(df)){
  set.seed(df$seed[r])
  x<-createData(sample.size=df$sample.size[r],gval=df$Gval[r],hval=df$Hval[r],contProp=df$contProp[r],contSpread=df$contSpread[r],error=as.character(df$error[r]))
  
  nsamp<-1000
  #Get the posterior

    BBmean.posterior<-BBmean(x,nsamp)
    BBmean.ci.low <- quantile(BBmean.posterior,probs=.025)
    BBmean.ci.high <- quantile(BBmean.posterior,probs=.975)


    BBtrim.posterior<-BBtrim(x,nsamp)
    BBtrim.ci.low <- quantile(BBtrim.posterior,probs=.025)
    BBtrim.ci.high <- quantile(BBtrim.posterior,probs=.975)

 
    normal.posterior<-normalModel(x,nsamp)
    normal.ci.low <- quantile(normal.posterior,probs=.025)
    normal.ci.high <- quantile(normal.posterior,probs=.975)

    cont.posterior<-CN.mean(x,nsamp)
    cont.ci.low <- cont.posterior$credible.interval[[1]]
    cont.ci.high <- cont.posterior$credible.interval[[2]]

  
  df$BBmeanWidth[r]<-abs(BBmean.ci.high-BBmean.ci.low)
  df$BBtrimWidth[r]<-abs(BBtrim.ci.high-BBtrim.ci.low)
  df$normalWidth[r]<-abs(normal.ci.high-normal.ci.low)
  df$contWidth[r]<-abs(cont.ci.high-cont.ci.low)
  df$no.contaminants[r]<-cont.posterior$number.contaminants
  remove(x)
  
  set.seed(df$seeds[r])
  contaminates<-rnorm(as.numeric(df$sample.size[r])*as.numeric(df$contProp[r]),0,df$contSpread[r])
  if (sum(contaminates>3)==length(contaminates) || sum(contaminates<-3)==length(contaminates)){
    df$contSkew[r]<-'YES'
  } else{
    df$contSkew[r]<-'NO'
  }
  
  if (BBmean.ci.low<0 & BBmean.ci.high>0){
    df$BBmeanCorrect[r]<-'YES'
  }else{
    df$BBmeanCorrect[r]<-'NO'
  }
  if (BBtrim.ci.low<0 & BBtrim.ci.high>0){
    df$BBtrimCorrect[r]<-'YES'
  }else{
    df$BBtrimCorrect[r]<-'NO'
  }
  if (normal.ci.low<0 & normal.ci.high>0){
    df$normalCorrect[r]<-'YES'
  }else{
    df$normalCorrect[r]<-'NO'
  }
  
  if (cont.ci.low<0 & cont.ci.high>0){
    df$contCorrect[r]<-'YES'
  }else{
    df$contCorrect[r]<-'NO'
  }
  
  print(r/nrow(df))
}



save(df,file="Simulations2_withContEst.Rdata")


colnames(df) <- c("Iteration" ,    "seeds"     , "contSkew",   "BBmean.Correct" ,"BBtrim.Correct", "normal.Correct" ,"cont.Correct" ,  "BBmean.Width"  , "BBtrim.Width" ,
                  "normal.Width"  , "cont.Width"  ,  "no.contaminants", "Gval"     ,     "Hval"       ,  "contProp"      ,"contSpread"   , "sample.size"  , "error"   )
df_long<-df %>% 
  gather(v, value, BBmean.Correct:cont.Width) %>% 
  separate(v, c("model", "col")) %>% 
  arrange(seeds) %>% 
  spread(col, value) %>%
ungroup()  
df_long$sample.size<-as.factor(df_long$sample.size)
df_long$error<-as.factor(df_long$error)
df_long$model<-as.factor(df_long$model)
df_long$Width<-as.numeric(df_long$Width)
df_long$contProp<-as.factor(df_long$contProp)

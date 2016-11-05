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

###################ERROR MODEL:##############################################
###################################################################################
                                #ACTUAL SCRIPT
####################################################################################
#df<-expand.grid(Iteration=1:100, correct='AAA',width=-1,Gval=c(0,0.5),Hval=c(0,0.5),contProp=c(0,.1,.2,.3,.4,.5),contSpread=c(2,5,10),sample.size=c(10,25,50,100,250),error=c('normal','normal.bias','unif','unif.bias','dp'),model=c('percentile_t','t_test'))
#df<-expand.grid(Iteration=1:100, correct='AAA',width=-1,Gval=0,Hval=0,contProp=c(0,.1,.2,.3,.4,.5),contSpread=2,sample.size=c(25,50,100,250),error=c('normal.bias','normal','dp'),model=c('BB mean','BB trim','normal','cont','percentile_t','t_test'))
df<-expand.grid(Iteration=1:100,seeds=-1,contSkew='AAA',BBmeanCorrect='AAA',BBtrimCorrect='AAA',normalCorrect='AAA',contCorrect='AAA',tCorrect='AAA',tWidth=-1,BBmeanWidth=-1,BBtrimWidth=-1,normalWidth=-1,contWidth=-1,BBmeanModes=-1,BBtrimModes=-1,normalModes=-1,contModes=-1,tModes=-1,no.contaminants=-1,Gval=0,Hval=0,contProp=c(0,.1,.2,.3,.4,.5),contSpread=10,sample.size=c(20,50,100,250),error=c('normal.bias','normal'))
rand.seed<-sample.int(size=nrow(df),n=100000,replace=FALSE)
df$BBmeanCorrect<-as.character(df$BBmeanCorrect)
df$BBtrimCorrect<-as.character(df$BBtrimCorrect)
df$normalCorrect<-as.character(df$normalCorrect)
df$contSkew<-as.character(df$contSkew)
df$contCorrect<-as.character(df$contCorrect)
df$tCorrect<-as.character(df$tCorrect)
df$seeds<-rand.seed
for (r in 1:nrow(df)){
  set.seed(df$seed[r])
  x<-createData(sample.size=df$sample.size[r],gval=df$Gval[r],hval=df$Hval[r],contProp=df$contProp[r],contSpread=df$contSpread[r],error=as.character(df$error[r]))
  
  nsamp<-1000
  #Get the posterior fo BBmean

    BBmean.posterior<-BBmean(x,nsamp)
    BBmean.HDI<-hdi(BBmean.posterior, allowSplit=TRUE,credMass=.95)
    if(is.null(nrow(BBmean.HDI))){ # IF Single Modal
      BBmean.ci.low <- BBmean.HDI[[1]]
      BBmean.ci.high <- BBmean.HDI[[2]]
      df$BBmeanWidth[r]<-abs(BBmean.ci.high-BBmean.ci.low)
      if (BBmean.ci.low<0 & BBmean.ci.high>0){
        df$BBmeanCorrect[r]<-'YES'
      }else{
        df$BBmeanCorrect[r]<-'NO'
      }
      df$BBmeanModes[r]<-1
    }else{ #if multi-modal posterior, calculate the accuracy and width of all sections
      BBmean.width.tmp<-0
      notcontained=TRUE
      for (i in 1:nrow(BBmean.HDI)){
        BBmean.ci.low.tmp <- BBmean.HDI[[i,1]]
        BBmean.ci.high.tmp <- BBmean.HDI[[i,2]]
        BBmean.width.tmp<-BBmean.width.tmp+abs(BBmean.ci.high.tmp-BBmean.ci.low.tmp)
        if(notcontained){
          if (BBmean.ci.low.tmp<0 & BBmean.ci.high.tmp>0){
            df$BBmeanCorrect[r]<-'YES'
            notcontained=FALSE
          }else{
            df$BBmeanCorrect[r]<-'NO'
          }
        }
      }
      df$BBmeanWidth[r]<-BBmean.width.tmp
      df$BBmeanModes[r]<-nrow(BBmean.HDI)
    }

##TRIMMED MODEL##    

    BBtrim.posterior<-BBtrim(x,nsamp)
    BBtrim.HDI<-hdi(BBtrim.posterior, allowSplit=TRUE,credMass=.95)
    if(is.null(nrow(BBtrim.HDI))){ # IF Single Modal
      BBtrim.ci.low <- BBtrim.HDI[[1]]
      BBtrim.ci.high <- BBtrim.HDI[[2]]
      df$BBtrimWidth[r]<-abs(BBtrim.ci.high-BBtrim.ci.low)
      if (BBtrim.ci.low<0 & BBtrim.ci.high>0){
        df$BBtrimCorrect[r]<-'YES'
      }else{
        df$BBtrimCorrect[r]<-'NO'
      }
      df$BBtrimModes[r]<-1
    }else{ #if multi-modal posterior, calculate the accuracy and width of all sections
      BBtrim.width.tmp<-0
      notcontained=TRUE
      for (i in 1:nrow(BBtrim.HDI)){
        BBtrim.ci.low.tmp <- BBtrim.HDI[[i,1]]
        BBtrim.ci.high.tmp <- BBtrim.HDI[[i,2]]
        BBtrim.width.tmp<-BBtrim.width.tmp+abs(BBtrim.ci.high.tmp-BBtrim.ci.low.tmp)
        if(notcontained){
          if (BBtrim.ci.low.tmp<0 & BBtrim.ci.high.tmp>0){
            df$BBtrimCorrect[r]<-'YES'
            notcontained=FALSE
          }else{
            df$BBtrimCorrect[r]<-'NO'
          }
        }
      }
      df$BBtrimWidth[r]<-BBtrim.width.tmp
      df$BBtrimModes[r]<-nrow(BBtrim.HDI)
    }
    
##NORMAL MODEL## 
    normal.posterior<-normalModel(x,nsamp)
    normal.HDI<-hdi(normal.posterior, allowSplit=TRUE,credMass=.95)
    if(is.null(nrow(normal.HDI))){ # IF Single Modal
      normal.ci.low <- normal.HDI[[1]]
      normal.ci.high <- normal.HDI[[2]]
      df$normalWidth[r]<-abs(normal.ci.high-normal.ci.low)
      if (normal.ci.low<0 & normal.ci.high>0){
        df$normalCorrect[r]<-'YES'
      }else{
        df$normalCorrect[r]<-'NO'
      }
      df$normalModes[r]<-1
    }else{ #if multi-modal posterior, calculate the accuracy and width of all sections
      normal.width.tmp<-0
      notcontained=TRUE
      for (i in 1:nrow(normal.HDI)){
        normal.ci.low.tmp <- normal.HDI[[i,1]]
        normal.ci.high.tmp <- normal.HDI[[i,2]]
        normal.width.tmp<-normal.width.tmp+abs(normal.ci.high.tmp-normal.ci.low.tmp)
        if(notcontained){
          if (normal.ci.low.tmp<0 & normal.ci.high.tmp>0){
            df$normalCorrect[r]<-'YES'
            notcontained=FALSE
          }else{
            df$normalCorrect[r]<-'NO'
          }
        }
      }
      df$normalWidth[r]<-normal.width.tmp
      df$normalModes[r]<-nrow(normal.HDI)
    }
    
    #contaminated model

    cont.out<-CN.mean(x,nsamp)
    cont.posterior<-cont.out$posterior
    cont.HDI<-hdi(cont.posterior, allowSplit=TRUE,credMass=.95)
    if(is.null(nrow(cont.HDI))){ # IF Single Modal
      cont.ci.low <- cont.HDI[[1]]
      cont.ci.high <- cont.HDI[[2]]
      df$contWidth[r]<-abs(cont.ci.high-cont.ci.low)
      if (cont.ci.low<0 & cont.ci.high>0){
        df$contCorrect[r]<-'YES'
      }else{
        df$contCorrect[r]<-'NO'
      }
      df$contModes[r]<-1
    }else{ #if multi-modal posterior, calculate the accuracy and width of all sections
      cont.width.tmp<-0
      notcontained=TRUE
      for (i in 1:nrow(cont.HDI)){
        cont.ci.low.tmp <- cont.HDI[[i,1]]
        cont.ci.high.tmp <- cont.HDI[[i,2]]
        cont.width.tmp<-cont.width.tmp+abs(cont.ci.high.tmp-cont.ci.low.tmp)
        if(notcontained){
          if (cont.ci.low.tmp<0 & cont.ci.high.tmp>0){
            df$contCorrect[r]<-'YES'
            notcontained=FALSE
          }else{
            df$contCorrect[r]<-'NO'
          }
        }
      }
      df$contWidth[r]<-cont.width.tmp
      df$contModes[r]<-nrow(cont.HDI)
    }
#t model
    
    t.posterior<-tModel(x,nsamp)
    t.HDI<-hdi(t.posterior, allowSplit=TRUE,credMass=.95)
    if(is.null(nrow(t.HDI))){ # IF Single Modal
      t.ci.low <- t.HDI[[1]]
      t.ci.high <- t.HDI[[2]]
      df$tWidth[r]<-abs(t.ci.high-t.ci.low)
      if (t.ci.low<0 & t.ci.high>0){
        df$tCorrect[r]<-'YES'
      }else{
        df$tCorrect[r]<-'NO'
      }
      df$tModes[r]<-1
    }else{ #if multi-modal posterior, calculate the accuracy and width of all sections
      t.width.tmp<-0
      notcontained=TRUE
      for (i in 1:nrow(t.HDI)){
        t.ci.low.tmp <- t.HDI[[i,1]]
        t.ci.high.tmp <- t.HDI[[i,2]]
        t.width.tmp<-t.width.tmp+abs(t.ci.high.tmp-t.ci.low.tmp)
        if(notcontained){
          if (t.ci.low.tmp<0 & t.ci.high.tmp>0){
            df$tCorrect[r]<-'YES'
            notcontained=FALSE
          }else{
            df$tCorrect[r]<-'NO'
          }
        }
      }
      df$tWidth[r]<-t.width.tmp
      df$tModes[r]<-nrow(t.HDI)
    }
  
  df$no.contaminants[r]<-cont.out$number.contaminants
  remove(x)
  
  set.seed(df$seeds[r])
  contaminates<-rnorm(as.numeric(df$sample.size[r])*as.numeric(df$contProp[r]),0,df$contSpread[r])
  if (sum(contaminates>3)==length(contaminates) || sum(contaminates<3)==length(contaminates)){
    df$contSkew[r]<-'YES'
  } else{
    df$contSkew[r]<-'NO'
  }
  
  print(r/nrow(df))
}



save(df,file=paste("Simulations2_withContEst9",10,".Rdata",sep=""))


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

# plot the simulation results 
load("~/Dropbox/Research/LaurenKennedy/RobustIntervals/lauren_simulations/SimulationsFinal.Rdata")
require(ggplot2)
require(lsr)
require(dplyr)
library(tidyr)

#Get it into long form (not wide form)
colnames(df) <- c("Iteration" ,    "seeds"     , "contSkew",   "BBmean.Correct" ,"BBtrim.Correct", "normal.Correct" ,"cont.Correct" ,  "BBmean.Width"  , "BBtrim.Width" ,
                  "normal.Width"  , "cont.Width"  ,  "no.contaminants", "Gval"     ,     "Hval"       ,  "contProp"      ,"contSpread"   , "sample.size"  , "error"   )
df_long<-df %>% 
  select(-no.contaminants) %>% #remove the column containing the number of contaminates the contaminated model removes 
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

### accuracy ###
df_long$group<-rep(NA,nrow(df_long))
df_long$group[df_long$Iteration<101]<-rep(1,19200)
df_long$group[df_long$Iteration>100 & df_long$Iteration<201]<-rep(2,19200)
df_long$group[df_long$Iteration>200 & df_long$Iteration<301]<-rep(3,19200)
df_long$group[df_long$Iteration>300 & df_long$Iteration<401]<-rep(4,19200)
df_long$group[df_long$Iteration>400 & df_long$Iteration<501]<-rep(5,19200)
df_long$group[df_long$Iteration>500 & df_long$Iteration<601]<-rep(6,19200)
df_long$group[df_long$Iteration>600 & df_long$Iteration<701]<-rep(7,19200)
df_long$group[df_long$Iteration>700 & df_long$Iteration<801]<-rep(8,19200)
df_long$group[df_long$Iteration>800 & df_long$Iteration<901]<-rep(9,19200)
df_long$group[df_long$Iteration>900 & df_long$Iteration<1001]<-rep(10,19200)
df_long$group<-as.factor(df_long$group)







  smalldata1 <- aggregate(Correct ~ error + model + sample.size + contProp + group, df_long, function(x){mean(x=="YES")})
smalldata<- do.call(data.frame, aggregate(. ~ error + model + sample.size + contProp , smalldata1, function(x){c(mean = mean(x), sd = sd(x))}))

smalldata$ci<-1.96*smalldata$Correct.sd/sqrt(10)

# rename, reorder and subset the models
levels(smalldata$model) <- c("BB mean", "BB trimmed",  "Contaminated","Normal")
smalldata <- smalldata[smalldata$model %in% c("BB mean","BB trimmed","Contaminated","Normal"), ]
smalldata$model <- droplevels(smalldata$model)
smalldata$model <- permuteLevels( smalldata$model, c(4,3,1,2))

# rename, reoder and subset the contaminant types
levels(smalldata$error) <- c("Biased normal","Unbiased normal")                                 
smalldata$error <- permuteLevels( smalldata$error, c(2,1)) 
smalldata <- smalldata[smalldata$error %in% c("Biased normal","Unbiased normal"),]
smalldata$error <- droplevels(smalldata$error)
levels(smalldata$error) <- c("Unbiased contaminants","Biased contaminants")

smalldata$marker <- NA
smalldata$marker[smalldata$model=="Normal"] <- 24
smalldata$marker[smalldata$model=="Contaminated"] <- 21
smalldata$marker[smalldata$model=="BB mean"] <- 17
smalldata$marker[smalldata$model=="BB trimmed"] <- 19

smalldata$contProp<-as.factor(smalldata$contProp)
# draw the picture
  p <- ggplot()

# p<-p+geom_errorbar(data=smalldata, 
#                  mapping=aes(x=sample.size, 
#                              group=model,
#                              ymin=Correct.mean-ci, 
#                              ymax=Correct.mean+ci), 
#                  width=.1) 

p <- p + geom_line(data=smalldata, 
                   mapping=aes(x=sample.size, 
                               y=Correct.mean, 
                               group=model))
p <- p + scale_shape_manual(values=c(24,21,17,19))
p <- p + geom_point(data=smalldata, 
                    mapping=aes(x=sample.size, 
                                y=Correct.mean,     
                                group=model,
                                shape=model),
                    size=3, fill="grey80")

p <- p + facet_grid(error ~ contProp)
p <- p + theme_bw()
p <- p + xlab("Sample Size") + ylab("Coverage Probability")
plot(p)

dev.print(pdf,file="simCoverCI.pdf", width=10,height=6)



### width ###
smalldata <- aggregate(Width ~ error + model + sample.size + contProp, df_long, mean)

# rename, reorder and subset the models
levels(smalldata$model) <- c("BB mean", "BB trimmed", "Normal", "Contaminated")
smalldata <- smalldata[smalldata$model %in% c("BB mean","BB trimmed", "Normal", "Contaminated"), ]
smalldata$model <- droplevels(smalldata$model)
smalldata$model <- permuteLevels( smalldata$model, c(3,4,1,2))

# rename, reoder and subset the contaminant types
levels(smalldata$error) <- c("Biased normal","Unbiased normal")                                 
smalldata$error <- permuteLevels( smalldata$error, c(2,1)) 
smalldata <- smalldata[smalldata$error %in% c("Biased normal","Unbiased normal"),]
smalldata$error <- droplevels(smalldata$error)
levels(smalldata$error) <- c("Unbiased contaminants","Biased contaminants")

# draw the picture
p <- ggplot()
p <- p + geom_line(data=smalldata, 
                   mapping=aes(x=sample.size, 
                               y=Width, 
                               group=model, 
                               show.legend=FALSE),
                   show.legend=FALSE)
p <- p + scale_shape_manual(values=c(24,21,17,19))
p <- p + geom_point(data=smalldata, 
                    mapping=aes(x=sample.size, 
                                y=Width,     
                                group=model, 
                                shape=model),
                    size=3, fill="grey80")
p <- p + facet_grid(error ~ contProp)
p <- p + theme_bw()
p <- p + xlab("Sample Size") + ylab("Average Interval Width")
plot(p)


dev.print(pdf,file="simWidth.pdf", width=10,height=6)


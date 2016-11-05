# plot the simulation results 
setwd("~/PhD/2016/Robust Write Up/Write Up/Reviewer Comments/HDI_Replication")

require(ggplot2)
require(lsr)
require(dplyr)
library(tidyr)

load("~/PhD/2016/Robust Write Up/Write Up/Reviewer Comments/HDI_Replication/Simulations2_withContEst10.Rdata")
df1<-df
load("~/PhD/2016/Robust Write Up/Write Up/Reviewer Comments/HDI_Replication/Simulations2_withContEst8.Rdata")
df2<-df
df2$Iteration<-df2$Iteration+100
load("~/PhD/2016/Robust Write Up/Write Up/Reviewer Comments/HDI_Replication/Simulations2_withContEst7.Rdata")
df3<-df
df3$Iteration<-df3$Iteration+200
load("~/PhD/2016/Robust Write Up/Write Up/Reviewer Comments/HDI_Replication/Simulations2_withContEst6.Rdata")
df4<-df
df4$Iteration<-df4$Iteration+300
load("~/PhD/2016/Robust Write Up/Write Up/Reviewer Comments/HDI_Replication/Simulations2_withContEst5.Rdata")
df5<-df
df5$Iteration<-df5$Iteration+400
load("~/PhD/2016/Robust Write Up/Write Up/Reviewer Comments/HDI_Replication/Simulations2_withContEst4.Rdata")
df6<-df
df6$Iteration<-df6$Iteration+500
load("~/PhD/2016/Robust Write Up/Write Up/Reviewer Comments/HDI_Replication/Simulations2_withContEst.Rdata")
df7<-df
df7$Iteration<-df7$Iteration+600
load("~/PhD/2016/Robust Write Up/Write Up/Reviewer Comments/HDI_Replication/Simulations2_withContEst3.Rdata")
df8<-df
df8$Iteration<-df8$Iteration+700
load("~/PhD/2016/Robust Write Up/Write Up/Reviewer Comments/HDI_Replication/Simulations2_withContEst2.Rdata")
df9<-df
df9$Iteration<-df9$Iteration+800
load("~/PhD/2016/Robust Write Up/Write Up/Reviewer Comments/HDI_Replication/Simulations2_withContEst9.Rdata")
df10<-df
df10$Iteration<-df9$Iteration+800
df<-rbind(df1,df2,df3,df4,df5,df6,df7,df8,df9,df10)




#Get it into long form (not wide form)
colnames(df) <- c("Iteration" ,    "seeds"     , "contSkew",   "BBmean.Correct" ,"BBtrim.Correct", "normal.Correct" ,"cont.Correct" , "t.Correct","t.Width", "BBmean.Width"  , "BBtrim.Width" ,
                  "normal.Width"  , "cont.Width"  , "BBmean.HDIModes", "BBtrim.HDIModes","normal.HDIModes","cont.HDIModes","t.HDIModes", "no.contaminants", "Gval"     ,     "Hval"       ,  
                  "contProp"      ,"contSpread"   , "sample.size"  , "error"   )
df_long<-df %>% 
  select(-no.contaminants) %>% #remove the column containing the number of contaminates the contaminated model removes 
  gather(v, value, BBmean.Correct:t.HDIModes) %>% 
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



library(dplyr)

smalldata<-df_long %>%
  group_by(contProp,sample.size,error,model) %>%
  mutate(mean.Width=mean(Width), proportion.contained=sum(Correct=='YES')/n())%>%
  ungroup()

# rename, reorder and subset the models

smalldata$model <- permuteLevels( smalldata$model, c(4,3,5,1,2))
levels(smalldata$error) <- c("Biased normal","Unbiased normal")                                 
smalldata$error <- permuteLevels( smalldata$error, c(2,1)) 
smalldata <- smalldata[smalldata$error %in% c("Biased normal","Unbiased normal"),]
smalldata$error <- droplevels(smalldata$error)
levels(smalldata$error) <- c("Unbiased contaminants","Biased contaminants")
# rename, reoder and subset the contaminant types

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
                               y=proportion.contained, 
                               group=model))
p <- p + scale_shape_manual(values=c(21,24,22,19,17))
p <- p + geom_point(data=smalldata, 
                    mapping=aes(x=sample.size, 
                                y=proportion.contained,     
                                group=model,
                                shape=model),
                    size=3, fill="grey80")

p <- p + facet_grid(error ~ contProp)
p <- p + theme_bw()
p <- p + xlab("Sample Size") + ylab("Coverage Probability")
plot(p)

dev.print(pdf,file="simCoverCI.pdf", width=10,height=6)



#width#

# draw the picture
p <- ggplot()
p <- p + geom_line(data=smalldata, 
                   mapping=aes(x=sample.size, 
                               y=mean.Width, 
                               group=model, 
                               show.legend=FALSE),
                   show.legend=FALSE)
p <- p + scale_shape_manual(values=c(21,24,22,19,17))
p <- p + geom_point(data=smalldata, 
                    mapping=aes(x=sample.size, 
                                y=mean.Width,     
                                group=model, 
                                shape=model),
                    size=3, fill="grey80")
p <- p + facet_grid(error ~ contProp)
p <- p + theme_bw()
p <- p + xlab("Sample Size") + ylab("Average Interval Width")
plot(p)


dev.print(pdf,file="simWidth.pdf", width=10,height=6)


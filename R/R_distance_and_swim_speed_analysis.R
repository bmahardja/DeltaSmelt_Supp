# Purpose: Swim distance plots and analyses
# Author: Brian Mahardja
# Date: 2023-11-27

library(tidyverse)
library(lubridate)
library(viridis)


# Set working directory
root <- "C:/Users/bmahardja/Documents/GitHub/DeltaSmelt_Supp"
setwd(root)

data_root<-file.path(root,"data")
code_root <- file.path(root,"R")
output_root <- file.path(root,"output")

data_WY22<-read.csv(file=file.path(output_root,"WY2022_distance_from_release_site.csv")) %>%
  mutate(SampleDate=as.Date(SampleDate)) %>% mutate(DaysinceNov1=as.numeric(SampleDate-as.Date("2021-11-01"))) %>%
  mutate(DaysinceNov1_sqr=DaysinceNov1^2)

data_WY23<-read.csv(file=file.path(output_root,"WY2023_distance_from_release_site.csv"))%>%
  mutate(SampleDate=as.Date(SampleDate)) %>% mutate(DaysinceNov1=as.numeric(SampleDate-as.Date("2022-11-01"))) %>%
  mutate(DaysinceNov1_sqr=DaysinceNov1^2)

mean(data_WY22$swim_speed)
mean(data_WY23$swim_speed)


plot_date<- ggplot(data=data_WY23,aes(x=SampleDate,y=swim_speed,color=ReleaseEvent)) + geom_point() +
  labs(title="WY23 Swim speed")
plot_date


tiff(filename=file.path(output_root,"Figure_date_swimspeed_WY23.tiff"),
     type="cairo",
     units="in", 
     width=7, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_date
dev.off()


plot_day<- ggplot(data=data_WY23,aes(x=DaySinceRelease,y=swim_speed,color=Release.Method)) + geom_point()+
  labs(title="WY23 Swim speed")
plot_day

tiff(filename=file.path(output_root,"Figure_daysincerelease_swimspeed_WY23.tiff"),
     type="cairo",
     units="in", 
     width=7, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_day
dev.off()

plot_date_22<- ggplot(data=data_WY22,aes(x=SampleDate,y=swim_speed)) + geom_point() +
  labs(title="WY22 Swim speed")
plot_date_22


tiff(filename=file.path(output_root,"Figure_date_swimspeed_WY22.tiff"),
     type="cairo",
     units="in", 
     width=7, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_date_22
dev.off()


plot_day_22<- ggplot(data=data_WY22,aes(x=DaySinceRelease,y=swim_speed)) + geom_point()+
  labs(title="WY22 Swim speed")
plot_day_22

tiff(filename=file.path(output_root,"Figure_daysincerelease_swimspeed_WY22.tiff"),
     type="cairo",
     units="in", 
     width=7, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_day_22
dev.off()


plot_dist_day_23<- ggplot(data=data_WY23,aes(x=DaySinceRelease,y=dist_km,color=Release.Method)) + geom_point() +
  labs(title="WY23 distance in km")
plot_dist_day_23

tiff(filename=file.path(output_root,"Figure_daysincerelease_dist_WY23.tiff"),
     type="cairo",
     units="in", 
     width=7, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_dist_day_23
dev.off()

plot_dist_day_22<- ggplot(data=data_WY22,aes(x=DaySinceRelease,y=dist_km)) + geom_point() +
  labs(title="WY22 distance in km")
plot_dist_day_22

tiff(filename=file.path(output_root,"Figure_daysincerelease_dist_WY22.tiff"),
     type="cairo",
     units="in", 
     width=7, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_dist_day_22
dev.off()

plot_dist_date_23<- ggplot(data=data_WY23,aes(x=SampleDate,y=dist_km,color=Release.Method)) + geom_point() +
  labs(title="WY23 distance in km")
plot_dist_date_23

tiff(filename=file.path(output_root,"Figure_date_dist_WY23.tiff"),
     type="cairo",
     units="in", 
     width=7, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_dist_date_23
dev.off()

plot_dist_date_22<- ggplot(data=data_WY22,aes(x=SampleDate,y=dist_km)) + geom_point() +
  labs(title="WY22 distance in km")
plot_dist_date_22

tiff(filename=file.path(output_root,"Figure_date_dist_WY22.tiff"),
     type="cairo",
     units="in", 
     width=7, #10*1, 
     height=6, #22*1, 
     pointsize=5, #12, 
     res=300,
     compression="lzw")
plot_dist_date_22
dev.off()





data_WY22 <- data_WY22 %>% mutate(log_swim_speed=log(swim_speed))


plot_date<- ggplot(data=data_WY22,aes(x=SampleDate,y=log_swim_speed)) + geom_point()
plot_date

plot_day<- ggplot(data=data_WY22,aes(x=DaySinceRelease,y=log_swim_speed)) + geom_point()
plot_day

plot_day<- ggplot(data=data_WY22,aes(x=DaySinceRelease,y=swim_speed)) + geom_point()
plot_day

plot_dist<- ggplot(data=data_WY22,aes(x=DaySinceRelease,y=dist_km)) + geom_point()
plot_dist

test<-lm(data=data_WY22,log_swim_speed~DaySinceRelease)
summary(test)

test<-lm(data=data_WY22,log_swim_speed~SampleDate)
summary(test)

test2<-glm(data=data_WY22,swim_speed~DaySinceRelease+SampleDate,family="poisson")
summary(test2)
plot(test2)

test3<-glm(data=data_WY22,dist_km~DaySinceRelease,family=Gamma(link="log"))
summary(test3)
plot(test3)

data_WY22$prediction<-predict(test3,newdata=data_WY22,type="response")
predict(test3,newdata=data_WY22,type="response")


plot_dist<- ggplot(data=data_WY22,aes(x=DaySinceRelease,y=dist_km)) + geom_point()
plot_dist

plot_dist<- ggplot(data=data_WY22,aes(x=SampleDate,y=dist_km)) + geom_point()
plot_dist

test4<-glm(data=data_WY22,swim_speed~DaySinceRelease,family=Gamma(link="log"))
summary(test4)
predict(test4,newdata=data_WY22,type="response")
exp(predict(test4,newdata=data_WY22,type="response"))

data_WY22$prediction<-exp(predict(test4,newdata=data_WY22,type="response"))

ggplot(data=data_WY22,aes(x=prediction,y=swim_speed)) + geom_point()
ggplot(data=data_WY22,aes(x=swim_speed,y=prediction)) + geom_point()

test5<-glm(data=data_WY23,swim_speed~DaySinceRelease,family=Gamma(link="identity"))

summary(test5)

data_WY23$prediction<-exp(predict(test5,newdata=data_WY23,type="response"))

ggplot(data=data_WY23,aes(x=swim_speed,y=prediction)) + geom_point()

###
test4<-glm(data=data_WY22,swim_speed~DaySinceRelease,family=Gamma(link="identity"))
test5<-glm(data=data_WY23,swim_speed~DaySinceRelease,family=Gamma(link="identity"))

test_lm1<- lm(log(swim_speed) ~ DaySinceRelease,data=data_WY22)
test_lm2<- lm(log(swim_speed) ~ DaySinceRelease,data=data_WY23)

summary(test_lm1)
summary(test_lm2)

test_lm3<- lm(log(swim_speed) ~ DaysinceNov1_sqr+DaysinceNov1,data=data_WY22)
test_lm4<- lm(log(swim_speed) ~ DaysinceNov1_sqr+DaysinceNov1,data=data_WY23)

summary(test_lm3)
summary(test_lm4)

##
example<-data.frame(DaySinceRelease=(c(1:110)))
example <-example %>% mutate(WY23pred=exp(predict(test5,newdata=example,type="response")),
                             WY22pred=exp(predict(test4,newdata=example,type="response"))) %>%
  gather(Year,"predicted_swimspeed",2:3)


ggplot() + geom_line(data=example,aes(x=DaySinceRelease,y=predicted_swimspeed,color=Year)) +
  geom_point(data = data_WY23,aes(x=DaySinceRelease,y=swim_speed))+
  geom_point(data = data_WY22,aes(x=DaySinceRelease,y=swim_speed),color="purple")

example2<-data.frame(DaySinceRelease=(c(1:110)))
example2 <-example2 %>% mutate(WY23pred=exp(predict(test_lm1,newdata=example2,type="response")),
                             WY22pred=exp(predict(test_lm2,newdata=example2,type="response"))) %>%
  gather(lm_yr,"predicted_swimspeed",2:3)


ggplot() + geom_line(data=example2,aes(x=DaySinceRelease,y=predicted_swimspeed,color=lm_yr)) +
  geom_point(data = data_WY23,aes(x=DaySinceRelease,y=swim_speed))+
  geom_point(data = data_WY22,aes(x=DaySinceRelease,y=swim_speed),color="purple")

exp(predict(test_lm1, example2, interval="predict"))


example3<-data.frame(DaysinceNov1=(c(40:150)),DaysinceNov1_sqr=(c(40:150)^2))
example3 <-example3 %>% mutate(WY22pred=exp(predict(test_lm3,newdata=example3,type="response")),
                               WY23pred=exp(predict(test_lm4,newdata=example3,type="response"))) %>%
  gather(lm_yr,"predicted_swimspeed",3:4)

ggplot() + geom_line(data=example3,aes(x=DaysinceNov1,y=predicted_swimspeed,color=lm_yr)) +
  geom_point(data = data_WY23,aes(x=DaysinceNov1,y=swim_speed))+
  geom_point(data = data_WY22,aes(x=DaysinceNov1,y=swim_speed),color="purple")


## EC50

rm(list= ls())
#
library(tidyverse)
library(readxl)
library(writexl)
library(directlabels)
library(tidyverse)
library(mixmeta)
library(dosresmeta)
library(splines)
library(ggplot2)

## setting wd
getwd()
setwd("~/Desktop/PollutionThreshold/")

## importing data
ec50 <- as.data.frame(read.csv("EC50_b.csv", header=TRUE))
str(ec50)

copper.fert <- ec50 %>% filter(Pollutant == "Copper") %>% 
  filter(Coral.response.s..measured == "FS")
plot(Threshold.Value ~ Value_est, data = copper, color = Species.Codes)
ggplot(copper.fert, 
       aes(x = Threshold.Estimate,
           y = Value_est)) + 
  geom_boxplot() 

copper.ls <- ec50 %>% filter(Pollutant == "Copper") %>% 
  filter(Coral.response.s..measured == "LS")
ggplot(copper.ls, 
       aes(x = Threshold.Estimate,
           y = Value_est)) + 
  geom_boxplot() 

copper.se <- ec50 %>% filter(Pollutant == "Copper") %>% 
  filter(Coral.response.s..measured == "SE")
ggplot(copper.se, 
       aes(x = Threshold.Estimate,
           y = Value_est)) + 
  geom_boxplot() 
 
nickel.fert <- ec50 %>% filter(Pollutant == "Nickel") %>% 
  filter(Coral.response.s..measured == "FS")
ggplot(nickel.fert, 
       aes(x = Threshold.Estimate,
           y = Value_est)) + 
  geom_boxplot() 

diuron.mqy <- ec50 %>% filter(Pollutant == "Diuron") %>% 
  filter(Coral.response.s..measured == "MQY")
ggplot(diuron.mqy, 
       aes(x = Threshold.Estimate,
           y = Value_est)) + 
  geom_boxplot() 

diuron.eqy <- ec50 %>% filter(Pollutant == "Diuron") %>% 
  filter(Coral.response.s..measured == "EQY")
ggplot(diuron.eqy, 
       aes(x = Threshold.Estimate,
           y = Value_est)) + 
  geom_boxplot()


## 9.30.20

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

## setting wd
getwd()
setwd("~/Desktop/PollutionThreshold/")

## importing data
poll <- as.data.frame(read.csv("Master Pollution Data_9.23.csv", header=TRUE))
str(poll)

### changing things to numeric & character as needed
poll <- poll %>% #converting appropriate columns to numeric...
  mutate_at(vars(level.converted.to.ug.L, N.for.computing.average, 
                 exposure.duration.in.days, Standardized.response.level, 
                 Standardized.SE, Standardized.SD, N.for.computing.average.1,
                 Time.to.response..numeric, Duration.of.response..numeric), as.character) %>% 
  mutate_at(vars(level.converted.to.ug.L, N.for.computing.average, 
                 exposure.duration.in.days, Standardized.response.level, 
                 Standardized.SE, Standardized.SD, N.for.computing.average.1,
                 Time.to.response..numeric, Duration.of.response..numeric), as.numeric)
str(poll) #"double"-check!

### calculating the cumulative exposure concentration at mg*day/L
poll$CumulativeExposure <- as.numeric(poll$level.converted.to.ug.L*(poll$exposure.duration.in.days))

### creating new column with RefID & Experiment 
### & creating new standard SD col with 0 instead of NA
poll <- poll %>%
  mutate(RefIDExp = paste(RefID, Experiment.Number, sep = "_")) %>% 
  mutate(Standardized.SD.NAis0 = coalesce(Standardized.SD, 0))


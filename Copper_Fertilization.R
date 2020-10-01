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

## creating another dataset without controls
poll.noControl <- poll %>% filter(CONTROL=="0")

############################
## creating subset datasets
############################
#### COPPER #####
#### EARLY ##
copper.fert <- poll %>% filter(Response.type == "fertilization success") %>% 
  filter(Pollutant == "Copper")
##
copper.sett <- poll %>% filter(Response.type == "Larval settlement") %>% 
  filter(Pollutant == "Copper")
##
copper.larvalsurvival <- poll %>% filter(Response.type == "Larval survival") %>% 
  filter(Pollutant == "Copper")
##
copper.all.early <- rbind(copper.fert, copper.sett, copper.larvalsurvival)

#### ADULT ##
copper.colsurv <- poll %>% filter(Response.type == "colony survival") %>% 
  filter(Pollutant == "Copper")
##
copper.bleach <- poll %>% filter(Response.type == "Bleaching") %>% 
  filter(Pollutant == "Copper")
##
copper.adult.blchmort <- rbind(copper.colsurv, copper.bleach)
##
##
copper.chla <- poll %>% filter(Response.type == "chl-a concentration") %>% 
  filter(Pollutant == "Copper")
##
copper.chlc <- poll %>% filter(Response.type == "chl-c concentration") %>% 
  filter(Pollutant == "Copper")
##
copper.chl_ac <- rbind(copper.chla, copper.chlc)
##
##
copper.maxquantyield <- poll %>% filter(Response.type == "maximum quantum yield - Fv/Fm") %>% 
  filter(Pollutant == "Copper")
##
copper.effquantyield <- poll %>% filter(Response.type == "effective quantum yield - delta F/Fm") %>% 
  filter(Pollutant == "Copper")
##
copper.photeff <- rbind(copper.maxquantyield, copper.effquantyield)


#### ALL METALS #####
#### EARLY ##
metals.fert <- poll %>% filter(Response.type == "fertilization success") %>% 
  filter(Pollutant.Class == "Metal")
##
metals.larvsurviv <- poll %>% filter(Response.type == "Larval survival") %>% 
  filter(Pollutant.Class == "Metal")
##
metals.larvsett <- poll %>% filter(Response.type == "Larval settlement") %>% 
  filter(Pollutant.Class == "Metal")
##
metals.allearly <- rbind(metals.fert, metals.larvsurviv, metals.larvsett)


#############################################################
### Using dosresmeta to calculate Hedges' d and variance of d
#############################################################
## and getting rid of experiments that don't have std dev
copper.fert.noNAinStDev <- copper.fert[complete.cases(copper.fert$Standardized.SD),]
### need to fix this later but Copper-09 exp2 doesnt have any data just a control 
### so for now am deleting it
copper.fert.noNAinStDev <- copper.fert.noNAinStDev[-c(98),] 

## & calculating Hedge's D & effect sizes
covar_copper.fert <- by(copper.fert.noNAinStDev, copper.fert.noNAinStDev$RefIDExp, function(x) 
  covar.smd(Standardized.response.level, Standardized.SD, N.for.computing.average.1, 
            "smd", method="hedges", data = x))
copper.fert.noNAinStDev$smd <- unlist(lapply(covar_copper.fert, function(x) x$y))
copper.fert.noNAinStDev$vmd <- unlist(lapply(covar_copper.fert, function(x) x$v))
copper.fert.noNAinStDev
## removing controls for modeling
copper.fert2 <- as.data.frame(subset(copper.fert.noNAinStDev, CONTROL=="0"))

## plots
copper.fert.cum <- ggplot(copper.fert2, 
                          aes(x = CumulativeExposure,
                              y = smd,
                              color = Species,
                              ymin = smd-vmd,
                              ymax = smd+vmd)) + 
  geom_pointrange() + geom_smooth(method = lm) +
  ggtitle("Copper - Fertilization success rate") +
  labs(x = "Copper cumulative exposure (ug x day/L)",
       y = "Effect size (Hedges' d +/- s.d., dosresmeta calculation)",
       color = "Species") +
  geom_abline(intercept=0, slope=0) +
  theme_classic() + scale_x_log10()
copper.fert.cum
#
copper.fert.conc <- ggplot(copper.fert2, 
                           aes(x = level.converted.to.ug.L,
                               y = smd,
                               color = Species,
                               ymin = smd-vmd,
                               ymax = smd+vmd)) + 
  geom_pointrange() + geom_smooth(method = lm) +
  ggtitle("Copper - Fertilization success rate") +
  labs(x = "Copper concentration (ug/L)",
       y = "Effect size (Hedges' d +/- s.d., dosresmeta calculation)",
       color = "Species") +
  geom_abline(intercept=0, slope=0) +
  theme_classic() + scale_x_log10()
copper.fert.conc
#########
copper.sett.noNAinStDev <- copper.sett[complete.cases(copper.sett$Standardized.SD),]
covar_copper.sett <- by(copper.sett.noNAinStDev, copper.sett.noNAinStDev$RefIDExp, function(x) 
  covar.smd(Standardized.response.level, Standardized.SD, N.for.computing.average.1, 
            "smd", method="hedges", data = x))
copper.sett.noNAinStDev$smd <- unlist(lapply(covar_copper.sett, function(x) x$y))
copper.sett.noNAinStDev$vmd <- unlist(lapply(covar_copper.sett, function(x) x$v))
copper.sett.noNAinStDev
copper.sett.noNAinStDev2 <- as.data.frame(subset(copper.sett.noNAinStDev, CONTROL=="0"))
#
copper.sett.cum <- ggplot(copper.sett.noNAinStDev2, 
                          aes(x = CumulativeExposure,
                              y = smd,
                              color = RefIDExp,
                              ymin = smd-vmd,
                              ymax = smd+vmd)) + 
  geom_pointrange() + geom_smooth(method = lm) +
  ggtitle("Copper - Larval Settlement Success rate") +
  labs(x = "Copper cumulative exposure (ug x day/L)",
       y = "Effect size (Hedges' d +/- s.d., dosresmeta calculation)",
       color = "RefIDEXP") +
  geom_abline(intercept=0, slope=0) +
  theme_classic() + scale_x_log10()
copper.sett.cum

##
copper.larvalsurvival.noNAinStDev <- copper.larvalsurvival[complete.cases(copper.larvalsurvival$Standardized.SD),]
covar_copper.larvalsurvival <- by(copper.larvalsurvival.noNAinStDev, copper.larvalsurvival.noNAinStDev$RefIDExp, function(x) 
  covar.smd(Standardized.response.level, Standardized.SD, N.for.computing.average.1, 
            "smd", method="hedges", data = x))
copper.larvalsurvival.noNAinStDev$smd <- unlist(lapply(covar_copper.larvalsurvival, function(x) x$y))
copper.larvalsurvival.noNAinStDev$vmd <- unlist(lapply(covar_copper.larvalsurvival, function(x) x$v))
copper.larvalsurvival.noNAinStDev
copper.larvalsurvival.noNAinStDev2 <- as.data.frame(subset(copper.larvalsurvival.noNAinStDev, CONTROL=="0"))
#
copper.larvalsurvival.cum <- ggplot(copper.larvalsurvival.noNAinStDev2, 
                                    aes(x = CumulativeExposure,
                                        y = smd,
                                        color = Species,
                                        ymin = smd-vmd,
                                        ymax = smd+vmd)) + 
  geom_pointrange() + geom_smooth(method = lm) +
  ggtitle("Copper - Larval survival rate") +
  labs(x = "Copper cumulative exposure (ug x day/L)",
       y = "Effect size (Hedges' d +/- s.d., dosresmeta calculation)",
       color = "Species") +
  geom_abline(intercept=0, slope=0) +
  theme_classic() + scale_x_log10()
copper.larvalsurvival.cum  
#####



copper.all.early
copper.all.early.noNAinStDev <- copper.all.early[complete.cases(copper.all.early$Standardized.SD),]
covar_copper.all.early <- by(copper.all.early.noNAinStDev, copper.all.early.noNAinStDev$RefIDExp, function(x) 
  covar.smd(Standardized.response.level, Standardized.SD, N.for.computing.average.1, 
            "smd", method="hedges", data = x))
copper.all.early.noNAinStDev$smd <- unlist(lapply(covar_copper.all.early, function(x) x$y))
copper.all.early.noNAinStDev$vmd <- unlist(lapply(covar_copper.all.early, function(x) x$v))
copper.all.early.noNAinStDev
copper.all.early.noNAinStDev2 <- as.data.frame(subset(copper.all.early.noNAinStDev, CONTROL=="0"))
#
copper.all.early.cum <- ggplot(copper.all.early.noNAinStDev2, 
                               aes(x = CumulativeExposure,
                                   y = smd,
                                   color = Species,
                                   ymin = smd-vmd,
                                   ymax = smd+vmd)) + 
  geom_pointrange() + geom_smooth(method = lm) +
  ggtitle("Copper - Fertilization success rate") +
  labs(x = "Copper cumulative exposure (ug x day/L)",
       y = "Effect size (Hedges' d +/- s.d., dosresmeta calculation)",
       color = "Species") +
  geom_abline(intercept=0, slope=0) +
  theme_classic() + scale_x_log10()
copper.all.early.cum




### model fitting 
plot(smd ~ CumulativeExposure, data = copper.fert2)
mod_copper.fert.1 <- mixmeta(smd ~ CumulativeExposure, S = vmd,
                             random =  ~ 1 | RefIDExp,
                             data = copper.fert2, method = "ml")
summary(mod_copper.fert.1)
hist(copper.fert2$smd)
hist(copper.fert2$vmd)
hist(copper.fert2$CumulativeExposure)

# LINEAR FIXED AND RANDOM EFFECTS NOT ACCOUNTING FOR WITHIN-STUDY CORRELATIONS
mod_copper.fert.a <- mixmeta(smd ~ CumulativeExposure, S = vmd,
                             random =  ~ CumulativeExposure | RefIDExp,
                             data = copper.fert2, method = "ml")
mod_copper.fert.a
summary(mod_copper.fert.a)
# LINEAR FIXED AND RANDOM, NESTED EFFECTS NOT ACCOUNTING FOR WITHIN-STUDY CORRELATIONS
mod_copper.fert.b <- mixmeta(smd ~ CumulativeExposure, S = vmd,
                             random =  ~ CumulativeExposure | RefID/RefIDExp,
                             data = copper.fert2, method = "ml")
mod_copper.fert.b
summary(mod_copper.fert.a)
# LINEAR FIXED AND RANDOM, NESTED EFFECTS NOT ACCOUNTING FOR WITHIN-STUDY CORRELATIONS
mod_copper.fert.c <- mixmeta(smd ~ CumulativeExposure, S = vmd,
                             random =  ~ CumulativeExposure | RefID/RefIDExp,
                             data = copper.fert2, method = "ml")
summary(mod_copper.fert.a)


# LINEAR FIXED AND RANDOM EFFECTS ACCOUNTING FOR WITHIN-COMPARISON CORRELATIONS
#Code for covariance matrix to include for addSlist for Comparison random effect
newlist_FERT <- list(NA)
for (i in seq(1,length(covar_copper.fert))) {
  newlist_FERT[i] <- list(covar_copper.fert[[i]]$S)
}
## using dataset without controls
unique(copper.fert2$RefIDExp)
mod_copper.fert.d <- mixmeta(smd ~ CumulativeExposure,
                             random = ~ 1 | RefIDExp,
                             data=copper.fert2, method="ml", #using dataset without controls
                             control=list(addSlist=newlist_FERT))




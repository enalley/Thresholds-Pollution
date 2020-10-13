### analyzing pollution data
## 9.23.20

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
library(tidyverse)
library(mixmeta)
library(dosresmeta)

setwd("~/Desktop/PollutionThreshold/")

pollute <- as.data.frame(read.csv("Master Pollution Data_10.08.csv", header=TRUE))
head(pollute)

### changing things to numeric & character as needed
pollute <- pollute %>% #converting appropriate columns to numeric...
  mutate_at(vars(level.converted.to.ug.L, N.for.computing.average, 
                 exposure.duration.in.days, Standardized.response.level, 
                 Standardized.SE, Standardized.SD, N.for.computing.average.1,
                 Time.to.response..numeric, Duration.of.response..numeric), as.character) %>% 
  mutate_at(vars(level.converted.to.ug.L, N.for.computing.average, 
                 exposure.duration.in.days, Standardized.response.level, 
                 Standardized.SE, Standardized.SD, N.for.computing.average.1,
                 Time.to.response..numeric, Duration.of.response..numeric), as.numeric)
#head(pollute) #"double"-check!

### creating new column with RefID & Experiment 
pollute <- pollute %>%
  mutate(RefIDExp = paste(RefID, Experiment.Number, sep = "_")) %>% 
  mutate(Gsp = paste(Current.genus.name, Current.species.name, sep = " ")) %>% 
  mutate(log10_ug_L = log10(level.converted.to.ug.L+1)) %>% 
  mutate(CumExp = level.converted.to.ug.L*exposure.duration.in.days) %>% # calculating the cumulative exposure concentration at mg*day/L
  select(RefID, RefIDExp, CONTROL, Pollutant.Class, Pollutant, Response.type, Current.genus.name, Gsp, Coral.age.class, level.converted.to.ug.L, log10_ug_L, CumExp, exposure.duration.in.days, Standardized.response.level, Standardized.SD, N.for.computing.average.1)

### DIURON-P03 needs to be removed
pollute <- pollute[!(pollute$RefID == "Diuron-P03"),]

## creating subsets
diuron.fert <- pollute %>% filter(Response.type == "fertilization success") %>% filter(Pollutant == "Diuron")
diuron.mqy <- pollute %>% filter(Response.type == "maximum quantum yield - Fv/Fm") %>% filter(Pollutant == "Diuron")
diuron.eqy <- pollute %>% filter(Response.type == "effective quantum yield - delta F/Fm") %>% filter(Pollutant == "Diuron")

## DIURON FERT
covar_diuron.fert <- by(diuron.fert, diuron.fert$RefIDExp, function(x) 
  covar.smd(Standardized.response.level, Standardized.SD, N.for.computing.average.1, 
            "smd", method="hedges", data = x))
diuron.fert$smd <- unlist(lapply(covar_diuron.fert, function(x) x$y))
diuron.fert$vmd <- unlist(lapply(covar_diuron.fert, function(x) x$v))

## DIURON MQY
covar_diuron.mqy <- by(diuron.mqy, diuron.mqy$RefIDExp, function(x) 
  covar.smd(Standardized.response.level, Standardized.SD, N.for.computing.average.1, 
            "smd", method="hedges", data = x))
diuron.mqy$smd <- unlist(lapply(covar_diuron.mqy, function(x) x$y))
diuron.mqy$vmd <- unlist(lapply(covar_diuron.mqy, function(x) x$v))

## DIURON EQY
covar_diuron.eqy <- by(diuron.eqy, diuron.eqy$RefIDExp, function(x) 
  covar.smd(Standardized.response.level, Standardized.SD, N.for.computing.average.1, 
            "smd", method="hedges", data = x))
diuron.eqy$smd <- unlist(lapply(covar_diuron.eqy, function(x) x$y))
diuron.eqy$vmd <- unlist(lapply(covar_diuron.eqy, function(x) x$v))

## removing controls for modeling
diuron.fert2 <- as.data.frame(subset(diuron.fert, CONTROL == "0"))
diuron.mqy2 <- as.data.frame(subset(diuron.mqy, CONTROL == "0"))
diuron.eqy2 <- as.data.frame(subset(diuron.eqy, CONTROL == "0"))

## PLOTS
#EXPOSURE CONCENTRATION vs. effect size, by study

### FERT
diuron.fert.conc <- ggplot(diuron.fert2, 
                           aes(x = level.converted.to.ug.L,
                               y = smd,
                               color = RefID,
                               ymin = smd-vmd,
                               ymax = smd+vmd)) + 
  geom_pointrange() + geom_smooth(method = lm) +
  ggtitle("Diuron - Fertilization success rate") +
  labs(x = "Diuron concentration (ug/L)",
       y = "Effect size (Hedges' d +/- s.d.)",
       color = "Study") +
  geom_abline(intercept=0, slope=0) +
  theme_classic() 
diuron.fert.conc
#
diuron.fert.conc2 <- ggplot(diuron.fert2, 
                            aes(x = level.converted.to.ug.L,
                                y = smd,
                                color = Gsp,
                                ymin = smd-vmd,
                                ymax = smd+vmd)) + 
  geom_pointrange() + geom_smooth(method = lm) +
  ggtitle("Diuron - Fertilization success rate") +
  labs(x = "Diuron concentration (ug/L)",
       y = "Effect size (Hedges' d +/- s.d.)",
       color = "Species") +
  geom_abline(intercept=0, slope=0) +
  theme_classic() + scale_x_log10()
diuron.fert.conc2

### MQY
diuron.mqy.conc <- ggplot(diuron.mqy2, 
                           aes(x = log10_ug_L,
                               y = smd,
                               color = RefID,
                               ymin = smd-vmd,
                               ymax = smd+vmd)) + 
  geom_pointrange() + geom_smooth(method = lm) +
  ggtitle("Diuron - MQY") +
  labs(x = "Diuron concentration (ug/L)",
       y = "Effect size (Hedges' d +/- s.d.)",
       color = "Study") +
  geom_abline(intercept=0, slope=0) +
  theme_classic() 
diuron.mqy.conc
#
diuron.mqy.conc2 <- ggplot(diuron.mqy2, 
                            aes(x = log10_ug_L,
                                y = smd,
                                color = Gsp,
                                ymin = smd-vmd,
                                ymax = smd+vmd)) + 
  geom_pointrange() + geom_smooth(method = lm) +
  ggtitle("Diuron - MQY") +
  labs(x = "Diuron concentration (ug/L)",
       y = "Effect size (Hedges' d +/- s.d.)",
       color = "Species") +
  geom_abline(intercept=0, slope=0) +
  theme_classic()
diuron.mqy.conc2

### EQY 
diuron.eqy.conc <- ggplot(diuron.eqy2, 
                          aes(x = log10_ug_L,
                              y = smd,
                              color = RefID,
                              ymin = smd-vmd,
                              ymax = smd+vmd)) + 
  geom_pointrange() + geom_smooth(method = lm) +
  ggtitle("Diuron - EQY") +
  labs(x = "Diuron concentration (ug/L)",
       y = "Effect size (Hedges' d +/- s.d.)",
       color = "Study") +
  geom_abline(intercept=0, slope=0) +
  theme_classic() 
diuron.eqy.conc
#
diuron.eqy.conc2 <- ggplot(diuron.eqy2, 
                           aes(x = log10_ug_L,
                               y = smd,
                               color = Gsp,
                               ymin = smd-vmd,
                               ymax = smd+vmd)) + 
  geom_pointrange() + geom_smooth(method = lm) +
  ggtitle("Diuron - MQY") +
  labs(x = "Diuron concentration (ug/L)",
       y = "Effect size (Hedges' d +/- s.d.)",
       color = "Species") +
  geom_abline(intercept=0, slope=0) +
  theme_classic()
diuron.eqy.conc2

### MODELS
# Organizational Note: 
#All models with 'a' at the end of its name have fixed slopes, random intercepts.
#All models with 'b' at the end of its name have random slopes and intercepts.
# LINEAR FIXED AND RANDOM EFFECTS NOT ACCOUNTING FOR WITHIN-STUDY CORRELATIONS
mod_diuron.fert.1a <- mixmeta(smd ~ level.converted.to.ug.L, S = vmd,
                              random =  ~ 1 | RefIDExp,
                              data = diuron.fert2, method = "ml")
summary(mod_diuron.fert.1a)
mod_diuron.fert.1b <- mixmeta(smd ~ level.converted.to.ug.L, S = vmd,
                              random =  ~ level.converted.to.ug.L | RefIDExp,
                              data = diuron.fert2, method = "ml")
summary(mod_diuron.fert.1b)

##
mod_diuron.mqy.1a <- mixmeta(smd ~ level.converted.to.ug.L, S = vmd,
                              random =  ~ 1 | RefIDExp,
                              data = diuron.mqy2, method = "ml")
summary(mod_diuron.mqy.1a)
mod_diuron.mqy.1b <- mixmeta(smd ~ level.converted.to.ug.L, S = vmd,
                              random =  ~ level.converted.to.ug.L | RefIDExp,
                              data = diuron.mqy2, method = "ml")
summary(mod_diuron.mqy.1b)

##
mod_diuron.eqy.1a <- mixmeta(smd ~ level.converted.to.ug.L, S = vmd,
                             random =  ~ 1 | RefIDExp,
                             data = diuron.eqy2, method = "ml")
summary(mod_diuron.eqy.1a)
mod_diuron.eqy.1b <- mixmeta(smd ~ level.converted.to.ug.L, S = vmd,
                             random =  ~ level.converted.to.ug.L | RefIDExp,
                             data = diuron.eqy2, method = "ml")
summary(mod_diuron.eqy.1b)

####
# LINEAR FIXED AND RANDOM EFFECTS ACCOUNTING FOR WITHIN-COMPARISON CORRELATIONS
#Code for covariance matrix to include for addSlist for RefIDExp random effect
newlist_FERT_diuron <- list(NA)
for (i in seq(1,length(covar_diuron.fert))) {
  newlist_FERT_diuron[i] <- list(covar_diuron.fert[[i]]$S)
}
mod_diuron.fert.3a <- mixmeta(smd ~ log10_ug_L,
                              random = ~ 1 | RefIDExp,
                              data=diuron.fert2, method="ml",
                              control=list(addSlist=newlist_FERT_diuron))
mod_diuron.fert.3b <- mixmeta(smd ~ log10_ug_L,
                              random = ~ log10_ug_L | RefIDExp,
                              data=diuron.fert2, method="ml",
                              control=list(addSlist=newlist_FERT_diuron))

##
# LINEAR FIXED AND RANDOM, NESTED EFFECTS ACCOUNTING FOR WITHIN-COMPARISON & WITHIN-STUDY CORRELATIONS
#Code for covariance matrix to include for addSlist for Ref/Comparison

### fertilization
newlist_FERT_diuron2 <- list(NA) #collects the list of covariance matrices for the block diag matrix by reference for the nested hierarchical ref/comparison
templist_FERT_diuron <- list(NA) #holds temp list of covariance matrices associated with each reference
templist_FERT_diuron2 <- list(NA)
reflista <- diuron.fert %>% distinct(RefID,RefIDExp)
reflist <- reflista[,1]
for (i in seq(1,length(unique(reflista$RefID))))  {
  #pull the elements from covar_copper.fert_b that are all from the same reference [i]
  templist_FERT_diuron[i] <-list(covar_diuron.fert[reflist==unique(reflista$RefID)[i]])
  for (j in seq(1,length(templist_FERT_diuron[[i]]))) {
    #for each comparison in the reference, pull out the covar matrices (element $S) and put in into templist_FERT_Cu2
    templist_FERT_diuron2[j] <- list(templist_FERT_diuron[[i]][[j]]$S)
  }
  #turn list of covars from all comparison in one reference into block diag matrix
  newlist_FERT_diuron2[i] <- list(bdiagMat(templist_FERT_diuron2))
  templist_FERT_diuron2 <- list(NA)
}
mod_diuron.fert.4a <- mixmeta(smd ~ level.converted.to.ug.L,
                              random = ~ 1 | RefID/RefIDExp,
                              data = diuron.fert2, method="ml",
                              control=list(addSlist=newlist_FERT_diuron2))
mod_diuron.fert.4b <- mixmeta(smd ~ level.converted.to.ug.L,
                              random = ~ level.converted.to.ug.L | RefID/RefIDExp,
                              data=diuron.fert2, method="ml",
                              control=list(addSlist=newlist_FERT_diuron2))


# Model Comparisons
AIC(mod_diuron.fert.1a, mod_diuron.fert.3a, 
    mod_diuron.fert.1b, mod_diuron.fert.3b) 
## 3a is the best -- the rest didn't work

# Residuals vs. Fitted Plots and Normal Q-Q Plot
op <- par(mfrow = c(1,2), mar = c(2,2,4,1))
resid_3a <- resid(mod_diuron.fert.3a)
fitted_3a <- fitted(mod_diuron.fert.3a)
plot(fitted_3a, resid_3a, main = "smd ~ log10_ug_L, S = vmd,
     random =  ~ log10_ug_L | RefID/RefIDExp")
abline(0,0)
qqnorm(resid_3a)
qqline(resid_3a)
par(op)

# PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS
pred_mod_diuron.fert.3a <- predict(mod_diuron.fert.3a, newdata=diuron.fert2, ci=TRUE)
diuron.fert2_CI <- cbind(diuron.fert2, pred_mod_diuron.fert.3a)
#head(copper.fert2_CI)
min_conc_diuron <- diuron.fert2_CI %>% 
  mutate(overlap0 = 0 >= ci.lb & 0 <= ci.ub) %>% 
  filter(overlap0==FALSE) %>% 
  filter(ci.ub<0) %>% 
  summarize(min_conc_diuron=min(level.converted.to.ug.L))
min_conc_diuron2 <- as.numeric(min_conc_diuron)

ggplot(diuron.fert2_CI, 
       aes(x = level.converted.to.ug.L,
           y = smd,
           color = RefID,
           ymin = smd-vmd,
           ymax = smd+vmd)) + 
  geom_pointrange() +
  labs(x = expression("Diuron exposure concentration ("*mu*"g/L)"),
       y = expression("Effect size (Hedges'"~italic(d)~"+/- variance of"~italic(d)~")"),
       color = "Study") +
  geom_abline(intercept=0, slope=0) +
  geom_vline(xintercept=min_conc_diuron2, linetype="dashed", color = "red") +
  geom_line(aes(x = level.converted.to.ug.L, y = fit), inherit.aes=FALSE) +
  geom_ribbon(aes(x = level.converted.to.ug.L, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  theme_classic() +
  scale_x_log10(limits = c(1,1000), breaks=c(1,10,100,1000,min_conc_diuron2), 
                label=c("1","","100","1000",round(min_conc_diuron2,digits=1)))


### mqy
newlist_MQY_diuron2 <- list(NA) #collects the list of covariance matrices for the block diag matrix by reference for the nested hierarchical ref/comparison
templist_MQY_diuron <- list(NA) #holds temp list of covariance matrices associated with each reference
templist_MQY_diuron2 <- list(NA)
reflista <- diuron.mqy %>% distinct(RefID,RefIDExp)
reflist <- reflista[,1]
for (i in seq(1,length(unique(reflista$RefID))))  {
  #pull the elements from covar_copper.fert_b that are all from the same reference [i]
  templist_MQY_diuron[i] <-list(covar_diuron.mqy[reflist==unique(reflista$RefID)[i]])
  for (j in seq(1,length(templist_MQY_diuron[[i]]))) {
    #for each comparison in the reference, pull out the covar matrices (element $S) and put in into templist_FERT_Cu2
    templist_MQY_diuron2[j] <- list(templist_MQY_diuron[[i]][[j]]$S)
  }
  #turn list of covars from all comparison in one reference into block diag matrix
  newlist_MQY_diuron2[i] <- list(bdiagMat(templist_MQY_diuron2))
  templist_MQY_diuron2 <- list(NA)
}
mod_diuron.mqy.4a <- mixmeta(smd ~ level.converted.to.ug.L,
                              random = ~ 1 | RefID/RefIDExp,
                              data = diuron.mqy2, method="ml",
                              control=list(addSlist=newlist_MQY_diuron2))
mod_diuron.mqy.4b <- mixmeta(smd ~ level.converted.to.ug.L,
                              random = ~ level.converted.to.ug.L | RefID/RefIDExp,
                              data=diuron.mqy2, method="ml",
                              control=list(addSlist=newlist_MQY_diuron2))
## not working!! 


### eqy
newlist_EQY_diuron2 <- list(NA) #collects the list of covariance matrices for the block diag matrix by reference for the nested hierarchical ref/comparison
templist_EQY_diuron <- list(NA) #holds temp list of covariance matrices associated with each reference
templist_EQY_diuron2 <- list(NA)
reflista <- diuron.eqy %>% distinct(RefID,RefIDExp)
reflist <- reflista[,1]
for (i in seq(1,length(unique(reflista$RefID))))  {
  #pull the elements from covar_copper.fert_b that are all from the same reference [i]
  templist_EQY_diuron[i] <-list(covar_diuron.eqy[reflist==unique(reflista$RefID)[i]])
  for (j in seq(1,length(templist_EQY_diuron[[i]]))) {
    #for each comparison in the reference, pull out the covar matrices (element $S) and put in into templist_FERT_Cu2
    templist_EQY_diuron2[j] <- list(templist_EQY_diuron[[i]][[j]]$S)
  }
  #turn list of covars from all comparison in one reference into block diag matrix
  newlist_EQY_diuron2[i] <- list(bdiagMat(templist_EQY_diuron2))
  templist_EQY_diuron2 <- list(NA)
}
mod_diuron.eqy.4a <- mixmeta(smd ~ level.converted.to.ug.L,
                             random = ~ 1 | RefID/RefIDExp,
                             data = diuron.eqy2, method="ml",
                             control=list(addSlist=newlist_EQY_diuron2))
mod_diuron.eqy.4b <- mixmeta(smd ~ level.converted.to.ug.L,
                             random = ~ level.converted.to.ug.L | RefID/RefIDExp,
                             data=diuron.eqy2, method="ml",
                             control=list(addSlist=newlist_EQY_diuron2))
# Model Comparisons
AIC(mod_diuron.eqy.1a, mod_diuron.eqy.1b, 
    mod_diuron.eqy.4a, mod_diuron.eqy.4b) 
## 1b is the best -- the rest didn't work

# Residuals vs. Fitted Plots and Normal Q-Q Plot
op <- par(mfrow = c(1,2), mar = c(2,2,4,1))
resid_1b <- resid(mod_diuron.eqy.1b)
fitted_1b <- fitted(mod_diuron.eqy.1b)
plot(fitted_1b, resid_1b, main = "smd ~ log10_ug_L, S = vmd,
     random =  ~ log10_ug_L | RefID/RefIDExp")
abline(0,0)
qqnorm(resid_1b)
qqline(resid_1b)
par(op)
## should drop value >1000

# PREDICT THE EFFECT SIZE AND PLOT WITH CONFIDENCE INTERVALS
pred_mod_diuron.eqy.1b <- predict(mod_diuron.eqy.1b, newdata=diuron.eqy2, ci=TRUE)
diuron.eqy2_CI <- cbind(diuron.eqy2, pred_mod_diuron.eqy.1b)
#head(copper.fert2_CI)
min_conc_diuron <- diuron.eqy2_CI %>% 
  mutate(overlap0 = 0 >= ci.lb & 0 <= ci.ub) %>% 
  filter(overlap0==FALSE) %>% 
  filter(ci.ub<0) %>% 
  summarize(min_conc_diuron=min(level.converted.to.ug.L))
min_conc_diuron2 <- as.numeric(min_conc_diuron)

ggplot(diuron.eqy2_CI, 
       aes(x = level.converted.to.ug.L,
           y = smd,
           color = RefID,
           ymin = smd-vmd,
           ymax = smd+vmd)) + 
  geom_pointrange() +
  labs(x = expression("Diuron exposure concentration ("*mu*"g/L)"),
       y = expression("Effect size (Hedges'"~italic(d)~"+/- variance of"~italic(d)~")"),
       color = "Study") +
  geom_abline(intercept=0, slope=0) +
  geom_vline(xintercept=min_conc_diuron2, linetype="dashed", color = "red") +
  geom_line(aes(x = level.converted.to.ug.L, y = fit), inherit.aes=FALSE) +
  geom_ribbon(aes(x = level.converted.to.ug.L, y = smd,
                  ymin = ci.lb, ymax = ci.ub), 
              fill = "grey70", alpha = .15, 
              show.legend=FALSE, inherit.aes=FALSE) +
  theme_classic() +
  scale_x_log10(limits = c(1,1000), breaks=c(1,10,100,1000,min_conc_diuron2), 
                label=c("1","","100","1000",round(min_conc_diuron2,digits=1)))
### outlier throwing off plot


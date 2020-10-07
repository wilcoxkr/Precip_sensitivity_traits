### 
### Running SEMs to assess linkages among traits and Sspp
###

### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
### Created: June 14, 2020; Last modified: June 14, 2020

###
### Set up workstation
###
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\data\\")

library(tidyverse)
library(lavaan)
library(lavaanPlot)

### 
### Read in data
###
source("..//scripts//02_compile ppt sens_traits_cluster.R")
rm(list=setdiff(ls(),"full_df"))

### Alternately, read in file directly
full_df <- read.csv("Sspp and traits_full.csv")

###
### Define focal traits for tidiness
###
trait_vec_sem <- c(
  "LogSLA",
  "LogLfOsmPot",
  "LogLfThickness",
  "LogLDMC",
  "LogLfN",
  "SqrtStemSpecDens"
)

###
### Put data in wide form with only focal traits for SEM
###
sem_df_wide <- full_df %>%
  drop_na() %>%
  dplyr::select(Species:life_form, LogSLA, LogLfOsmPot, LogLfThickness, LogLDMC, LogLfN, SqrtStemSpecDens)  %>%
  mutate(ppt_slope_sqrt = sqrt(ppt_slope + abs(min(ppt_slope)) + .001) )

### Check correlation between LDMC and SLA and other variables because of potential collinearity issues
ggplot(full_df, aes(x=LogSLA, y=LogLDMC)) + geom_point() + theme_bw() +
  annotate("text", x=2.3, y=-0.3, label=paste0("r=",cor(full_df$LogLDMC, full_df$LogSLA)))
with(sem_df_wide, cor(LogSLA, LogLfOsmPot))
with(sem_df_wide, cor(LogLfThickness, LogLfOsmPot))
with(sem_df_wide, cor(LogSLA, LogLfOsmPot))
scatterplotMatrix(~ ppt_slope_sqrt + LogSLA + LogLfOsmPot + LogLfThickness + LogLDMC, data=filter(sem_df_wide, Species!="LIIN"),
                  regLine=T, smooth=F)

###
### Construct models for SEM
###

apriori_pa_model <- '
  ppt_slope_sqrt ~ LogLfOsmPot + LogLDMC + LogSLA + LogLfThickness
  LogLfOsmPot ~ LogLDMC
  LogSLA ~ LogLDMC + LogLfThickness

  LogLfThickness ~~ LogLDMC
 '

### Testing additional models based on reviewer comments
alternate1_pa_model <- '
  ppt_slope_sqrt ~ LogLfOsmPot + LogLDMC + LogSLA + LogLfThickness
  LogLDMC ~ LogLfOsmPot
  LogSLA ~ LogLDMC + LogLfThickness

  
  LogLfThickness ~~ LogLDMC
 '

alternate3_pa_model <- '
  ppt_slope_sqrt ~ LogLfOsmPot + LogLDMC + LogSLA
  LogLDMC ~ LogLfOsmPot
  LogSLA ~ LogLDMC
 '

alternate4_pa_model <- '
  ppt_slope_sqrt ~ LogLfOsmPot + LogLDMC + LogSLA + LogLfThickness
  LogLfOsmPot ~ LogLDMC
  LogSLA ~ LogLDMC + LogLfThickness

  
  LogLfThickness ~~ LogLDMC
 '

alternate5_pa_model <- '
  ppt_slope_sqrt ~ LogLDMC + LogSLA + LogLfThickness
  LogSLA ~ LogLDMC + LogLfThickness

  LogLfThickness ~~ LogLDMC
 '


###
### Fit models and print summaries
###

apriori_pa_model_fit <- sem(apriori_pa_model, data=sem_df_wide) # Best fit of data to model
summary(apriori_pa_model_fit, standardized=TRUE)
parameterEstimates(apriori_pa_model_fit)
AIC(apriori_pa_model_fit)

alternate1_pa_model_fit <- sem(alternate1_pa_model, data=sem_df_wide) 
summary(alternate1_pa_model_fit, standardized=TRUE)
parameterEstimates(alternate1_pa_model_fit)
AIC(alternate1_pa_model_fit)

alternate3_pa_model_fit <- sem(alternate3_pa_model, data=sem_df_wide) 
summary(alternate3_pa_model_fit, standardized=TRUE)
parameterEstimates(alternate3_pa_model_fit)
AIC(alternate3_pa_model_fit)

alternate4_pa_model_fit <- sem(alternate4_pa_model, data=sem_df_wide) 
summary(alternate4_pa_model_fit, standardized=TRUE)
parameterEstimates(alternate4_pa_model_fit)
AIC(alternate4_pa_model_fit)

alternate5_pa_model_fit <- sem(alternate5_pa_model, data=sem_df_wide) 
summary(alternate5_pa_model_fit, standardized=TRUE)
parameterEstimates(alternate5_pa_model_fit)
AIC(alternate5_pa_model_fit)



###
### Plotting
###

## Plotted in inkscape based on Std.all column from apriori_pa_model_fit summary



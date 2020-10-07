### Compile sensitivity (ppt) + fxn groups + trait dataset 
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last updated June 14, 2020

library(tidyverse)

#setwd("T:\\1-People\\Dana Lab\\TRAITS\\PUBLICATIONS\\precip_sensitivity\\")
#setwd("C:\\Users\\wilco\\Desktop\\from ARS comp\\precip_sensitivity\\")
#setwd("C:\\Users\\kwilcox4\\Dropbox\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # work computer
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # HP laptop

source("scripts\\01_calculate precipitation sensitivity with standardized covers.R")

fxn_groups <- read.csv("data//fxn_groups_cper.csv")

sp_with_traits <- read.csv("data\\trait_data\\Mean trait values summary table_GOLD_Tranformed.csv") %>%
  filter(Site=="CPER")

full_df <- mean_sens_rel %>%
  rename(ppt_slope=rel_slope_mean, ppt_slope_se=rel_slope_se_mean) %>%
  dplyr::select(-rownum) %>%
  left_join(fxn_groups, by=c("species"="sp_code")) %>%
  full_join(sp_with_traits, by=c("species"="Species")) %>%
  rename(cluster=cluster_bray_with_volume_transformed) %>%
  rename(Species=species) %>%
  dplyr::select(-Site, 
                -Functional.Group, 
                -Phs_area, 
                -LogLfWatPot, 
                -SqrtStomCond,
                -LogRtDiam,
                -LogRDMC,
                -LogSRL,
                -SqrtRtN)
  


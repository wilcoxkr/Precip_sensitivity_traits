### LTGI foliar cover prep

### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last modified: June. 19, 2020

library(tidyverse)
library(ggplot2)
library(ggthemes)

#setwd("C:\\Users\\wilco\\Dropbox\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\")
#setwd("C:\\Users\\kwilcox4\\Dropbox\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # work desktop
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # HP laptop

plots_2_exclude <- read.csv("data\\cover_data\\LTGI_BasCov_wFlags.csv") %>%
  dplyr::select(Year, Pasture, Plot, Flag) %>%
  filter(Flag != "")

sp_with_traits <- read.csv("data\\trait_data\\Mean trait values summary table_GOLD_Tranformed.csv") %>%
  filter(Site=="CPER")

sp_keepers_cper <- as.character(sp_with_traits$Species)

ltgi_pasture_key <- read.csv("data\\cover_data\\LTGI_pasture_key.csv")

ltgi_cover_full <-  read.csv("data\\cover_data\\LTGI_Foliar_abs_cover_Feb2016.csv") %>%
  dplyr::select(-BARE,
                -DUNG,
                -LICH,
                -LITT,
                -MOSS,
                -TOTAL) %>%
  rename(Year=SamplingYear) %>%
  left_join(plots_2_exclude, by=c("Year","Pasture","Plot")) %>%
  filter(is.na(Flag)) %>%
  dplyr::select(-Flag) %>%
  left_join(ltgi_pasture_key, by="Pasture") %>%
  mutate(trt_grp = ifelse(Treatment %in% c("UN","LG"), "ULG", "MHG")) # lumps ungrazed with light grazed & mod with heavy grazed

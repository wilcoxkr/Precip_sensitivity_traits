### Data cleaning for GZTX basal cover data
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last modified: Oct. 10, 2017

library(tidyverse)
#setwd("T:\\1-People\\Dana Lab\\TRAITS\\PUBLICATIONS\\precip_sensitivity\\")
#setwd("C:\\Users\\wilco\\Desktop\\from ARS comp\\precip_sensitivity\\")
#setwd("C:\\Users\\kwilcox4\\Dropbox\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # work computer
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\")

ID_2_remove <- c(43046, 44000, 93507)

gztx_cover_raw <- read.csv("data\\cover_data\\GZTX_BasCov_LTEAR_1992_2017.csv") %>%
  dplyr::select(ID, Pasture, SamplingYear, Treatment, Quadrat, Species, New_Mid.Point) %>%
  rename(Year=SamplingYear, abs_cover=New_Mid.Point) %>%
  mutate(Species=toupper(Species)) %>%
  filter(!Species %in% c("BARE","DUNG","FUNGI","LICH","LIT","LITT","MUSH","MUSHI","SCAT","OPPO")) %>%
  filter(Treatment %in% c("UNUN","GZGZ")) %>%
  filter(!ID %in% ID_2_remove) %>%
  mutate(Pasture=replace(Pasture, ID==59640, 19)) %>%
  mutate(Species=replace(Species, Species=="HECO", "STCO")) %>%
  mutate(Species=replace(Species, Species=="AGSM", "PASM"))





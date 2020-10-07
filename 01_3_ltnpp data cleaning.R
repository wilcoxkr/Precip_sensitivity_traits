### preparing/cleaning LTNPP datasets
### datasets included are MIDSLOPE & RIDGE (83-08, 14-16), ESA (83-08), SEC25 (83-08)
###
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last modified: Dec. 6, 2017

library(tidyverse)
library(ggplot2)
library(ggthemes)

#setwd("T:\\1-People\\Dana Lab\\TRAITS\\PUBLICATIONS\\precip_sensitivity\\")
#setwd("C:\\Users\\wilco\\Desktop\\from ARS comp\\precip_sensitivity\\")
# setwd("C:\\Users\\kwilcox4\\Dropbox\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # work computer
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\")


### NPP wasn't collected to species level in 2009-2012, so I've removed those years from species level data
### I think the best way to proceed is to remove these years for all calculations in this dataset
### Also, 2013 data do not exist so we're missing 2009-2013
### IMPORTANT NOTE: anpp.csv file has raw weights from subplots... all weights need to be multiplied by 4 to convert to g m-2

### Duplicate species in dataset include:
### "CADU" and "CAEL" and "CAHE" - all three are Carex duriuscula -- (in ltnpp) (using CAEL nomenclature)
### "PASM" and "AGSM" - both Pascopyrum smithii -- (in ltnpp)
### "ELEL" and "SIHY" - both Elymus elymoides -- (in ltnpp) (only SIHY in esa/sec25) (using SIHY nomenclature)
### "HECO" and "STCO" - Hesperostipa comata -- (in ltnpp) (only STCO in esa/sec25)(using STCO nomenclature)
### "CHVI" and "HEVI" - Heterotheca villosa -- (in ltnpp) (only CHVI in esa/sec25)
### "SAIB" and "SAKA" - Salsola iberica -- (in ltnpp) (both in esa/sec25)
### "SPCO" and "Spco" - Sphaealcea coccinea (both in esa/sec25)
### "THFI" and "THTR" - Thelesperma filifolium (both in esa/sec25)
### "PLPA  " and "PLPA" - both Plantago patagonica (both in esa/sec25)

### Identify species to remove (shrubs)
shrub_vec <- c("ARFR", "ATCA", "CELA", "CHNA", "EREF", "GUSA")

### MIDSLOPE, RIDGE and SWALE data
ltnpp_base_long <-  read.csv("data\\npp_data\\CPER_LTNPP_spNPP_83-16.csv") %>%
  filter(site_name %in% c("MIDSLOPE","RIDGE", "SWALE")) %>%
  filter(!Year %in% 2009:2012) %>%
  mutate(CAEL_merged = ifelse(is.na(CAEL),CAHE, CAEL)) %>%
  mutate(CADU_merged = ifelse(is.na(CADU),CAEL_merged, CADU)) %>%
  mutate(PASM_merged = ifelse(is.na(PASM),AGSM, PASM)) %>%
  mutate(ELEL_merged = ifelse(is.na(ELEL),SIHY, ELEL)) %>%
  mutate(HECO_merged = ifelse(is.na(HECO),STCO, HECO)) %>%
  mutate(HEVI_merged = ifelse(is.na(HEVI),CHVI, HEVI)) %>%
  mutate(SAIB_merged = ifelse(is.na(SAIB),SAKA, SAIB)) %>%
  dplyr::select(-CAEL,# combined CAEL and CAHE into one column
                -CAHE,
                -CAEL_merged,
                -CADU, # combined CADU and CAEL_merged into one column
                -PASM,# combined AGSM and PASM into one colum
                -AGSM,
                -ELEL,# combined ELEL and SIHY
                -SIHY,
                -HECO,# combined HECO and STCO
                -STCO,
                -HEVI,# combined HEVI and CHVI
                -CHVI,
                -SAIB,# combined SAIB and SAKA 
                -SAKA,
                -attribute_LKU,
                -OSD, # old standing dead -- previous year growth
                -CSAG, # cool season annual gram
                -CSPG, # cool season perennial gram
                -FORB, # forbs
                -BOBU, # bouteloua gracilis + bouteloua dacteloides
                -WSPG, # warm season perennial gram
                -SS # sub shrub
                ) %>%
  rename(CAEL=CADU_merged, PASM=PASM_merged, SIHY=ELEL_merged, STCO=HECO_merged,
         HEVI=HEVI_merged, SAIB=SAIB_merged, year=Year) %>%
  gather(key=species, value=anpp, -(year:plot)) %>%
#  filter(!species %in% shrub_vec) %>% ## Removes ARFR from dataset... this is a test to see if the presence of shrubs is artifically decreasing sensitivity (because shrubs are hard to measure)
  replace_na(list(anpp=0))

### ESA and SEC25 data
esa_sec25_df_long <- read.csv("data\\npp_data\\anpp.csv") %>%
  filter(!comments %in% c("Duplicate Bag #2","Duplicate_krw")) %>%
  dplyr::select(-(attribute_LKU:Link_ID),-anpp_ID) %>%
  filter(site_name %in% c("ESA", "SEC25")) %>%
  filter(!year %in% 2009:2012) %>%
  mutate(species=replace(species, species=="AGSM", "PASM")) %>%
  mutate(species=replace(species, species=="ELEL", "SIHY")) %>%
  mutate(species=replace(species, species=="HECO", "STCO")) %>%
  mutate(species=replace(species, species=="CHVI", "HEVI")) %>%
  mutate(species=replace(species, species=="SAKA", "SAIB")) %>%
  mutate(species=replace(species, species=="Spco", "SPCO")) %>%
  mutate(species=replace(species, species=="THTR", "THFI")) %>%
  mutate(species=replace(species, species=="PLPA  ", "PLPA")) %>%
  mutate(species=replace(species, species=="CADU", "CAEL")) %>%
  mutate(species=replace(species, species=="CAHE", "CAEL")) %>%
#  filter(!species %in% shrub_vec) %>% ## Removes ARFR from dataset... this is a test to see if the presence of shrubs is artifically decreasing sensitivity (because shrubs are hard to measure)
  spread(key=species, value=weight) %>%
  # mutate(CAEL_merged = ifelse(is.na(CAEL),CAHE, CAEL)) %>%
  # dplyr::select(-CAEL,# combined CAEL and CAHE into one column
  #               -CAHE
  #               ) %>%
  # rename(CAEL=CAEL_merged) %>%
  gather(key=species, value=anpp, -(year:plot)) %>%
  replace_na(list(anpp=0)) %>%
  mutate(anpp=anpp*4)

### Combine datasets
ltnpp_sp_npp <- rbind(esa_sec25_df_long, ltnpp_base_long)

### Calculate plot NPP
ltnpp_plot_npp <- ltnpp_sp_npp %>%
  group_by(year, site_name, transect, plot) %>%
  summarise(anpp=sum(anpp,na.rm=T))

### Calculate pasture NPP
ltnpp_pasture_npp <- ltnpp_plot_npp %>%
  group_by(year, site_name) %>%
  summarise(anpp=mean(anpp,na.rm=T))

### Just ESA without ARFR
esa_noarfr_df_long <- read.csv("data\\npp_data\\anpp.csv") %>%
  filter(!comments %in% c("Duplicate Bag #2","Duplicate_krw")) %>%
  dplyr::select(-(attribute_LKU:Link_ID),-anpp_ID) %>%
  filter(site_name %in% c("ESA")) %>%
  filter(!year %in% 2009:2012) %>%
  mutate(species=replace(species, species=="AGSM", "PASM")) %>%
  mutate(species=replace(species, species=="ELEL", "SIHY")) %>%
  mutate(species=replace(species, species=="HECO", "STCO")) %>%
  mutate(species=replace(species, species=="CHVI", "HEVI")) %>%
  mutate(species=replace(species, species=="SAKA", "SAIB")) %>%
  mutate(species=replace(species, species=="Spco", "SPCO")) %>%
  mutate(species=replace(species, species=="THTR", "THFI")) %>%
  mutate(species=replace(species, species=="PLPA  ", "PLPA")) %>%
  mutate(species=replace(species, species=="CADU", "CAEL")) %>%
  mutate(species=replace(species, species=="CAHE", "CAEL")) %>%
  filter(species != "ARFR") %>% ## Removes ARFR from dataset... this is a test to see if the presence of shrubs is artifically decreasing sensitivity (because shrubs are hard to measure)
  spread(key=species, value=weight) %>%
  # mutate(CAEL_merged = ifelse(is.na(CAEL),CAHE, CAEL)) %>%
  # dplyr::select(-CAEL,# combined CAEL and CAHE into one column
  #               -CAHE
  #               ) %>%
  # rename(CAEL=CAEL_merged) %>%
  gather(key=species, value=anpp, -(year:plot)) %>%
  replace_na(list(anpp=0)) %>%
  mutate(anpp=anpp*4)

### Calculate plot NPP
esa_noarfr_plot_npp <- esa_noarfr_df_long %>%
  group_by(year, site_name, transect, plot) %>%
  summarise(anpp=sum(anpp,na.rm=T))

### Calculate pasture NPP
esa_noarfr_pasture_npp <- esa_noarfr_plot_npp %>%
  group_by(year, site_name) %>%
  summarise(anpp=mean(anpp,na.rm=T))


rm(esa_sec25_df_long, ltnpp_base_long)




### Calculating precipitation senstivity for CPER (LTGI)
###   Binning ungrazed with light grazing and moderate with heavy grazing
###   For generating precipitation sensitivity data set

### Author: Kevin Wilcox (kevin.wilcox@uwyo.edu)
###
### Last modified: Aug. 28, 2019

library(tidyverse)
library(ggthemes)
library(car)

#setwd("T:\\1-People\\Dana Lab\\TRAITS\\PUBLICATIONS\\precip_sensitivity\\")
#setwd("C:\\Users\\wilco\\Desktop\\from ARS comp\\precip_sensitivity\\")
# setwd("C:\\Users\\kwilcox4\\Dropbox\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # work desktop
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

# Read in cper precip data and get water year precip ----------------------
ppt_water_year <- read.csv("data\\ppt_data\\CPE_dailyprecip_1939-2017Sept.csv") %>%
  mutate(ppt_mm = INCHES*25.4) %>%
  dplyr::select(-INCHES) %>%
  rename(year=YEAR,month=MONTH,day=DAY) %>%
  mutate(water_year = ifelse(month<=8, year, year+1)) %>%
  group_by(water_year) %>%
  summarise(water_year_ppt = sum(ppt_mm,na.rm=T)) %>%
  filter(water_year %in% c(1992:2017)) %>%
  rename(year=water_year)

ppt_annual <- read.csv("data\\ppt_data\\CPE_dailyprecip_1939-2017Sept.csv") %>%
  mutate(ppt_mm = INCHES*25.4) %>%
  dplyr::select(-INCHES) %>%
  rename(year=YEAR,month=MONTH,day=DAY) %>%
  group_by(year) %>%
  summarise(annual_ppt = sum(ppt_mm,na.rm=T))
#sd(ppt_annual$annual_ppt)

# Create total plot cover and cover by species datasets
ltgi_pasture_tot_cover_ppt <- ltgi_cover_full %>%
  mutate(total_cover = rowSums(ltgi_cover_full[5:ncol(ltgi_cover_full)-2])) %>%
  dplyr::select(Year, Pasture, Transect, Plot, total_cover) %>%
  group_by(Year, Pasture) %>%
  summarise(total_cover=mean(total_cover, na.rm=T)) %>%
  left_join(ppt_water_year, by=c("Year"="year"))

ltgi_pasture_sp_cover_ppt <- ltgi_cover_full %>%
  gather(key=species, value=abs_cover, -Year, -Pasture, -Transect, -Plot, -Treatment, -trt_grp) %>%
  group_by(Year, Pasture, species) %>%
  summarise(abs_cover=mean(abs_cover, na.rm=T)) %>%
  left_join(ppt_water_year, by=c("Year"="year")) %>%
  mutate(trt_grp = ifelse(Pasture %in% c("GS","23W"), "ULG", "MHG")) %>% # lumps ungrazed with light grazed & mod with heavy grazed
  ungroup() %>%
  group_by(trt_grp, species) %>%
  mutate(abs_cover_pchange = (abs_cover-mean(abs_cover))/mean(abs_cover))

ltgi_pasture_sp_cover_ppt$abs_cover_pchange[is.nan(ltgi_pasture_sp_cover_ppt$abs_cover_pchange)] <- 0  


# calculate ppt-cover slopes by species and relativize slopes by mean cover in each treatment group
trt_grp_vec <- unique(ltgi_pasture_sp_cover_ppt$trt_grp)

Sens_out_water_year <- {}

num_nonzeros <- function(x){length(x[x > 0])} # function for identifying non-0 entries

no_cover_sp <- data.frame(species=levels(factor(ltgi_pasture_sp_cover_ppt$species)),ltgi=1) %>%
  full_join(data.frame(species=sp_keepers_cper,trait=1), by="species") %>%
  filter(is.na(ltgi)) %>%
  dplyr::select(species)

for(j in 1:length(trt_grp_vec)){
  
  # remove species without cover
  sum_covers <- ltgi_pasture_sp_cover_ppt %>%
    filter(trt_grp==trt_grp_vec[j]) %>%
    group_by(species) %>%
    summarise(sum_cov=sum(abs_cover),
              no_entries=num_nonzeros(abs_cover))
  
  sp_exclude <- c(subset(sum_covers, sum_cov < 0.000000001 | no_entries <= 1)$species,no_cover_sp$species)
  sp_vec <- sp_keepers_cper[!sp_keepers_cper %in% sp_exclude]
  #rm(sum_covers, sp_exclude)
  
  for(i in 1:length(sp_vec)){
    cover_temp <- subset(ltgi_pasture_sp_cover_ppt, species==sp_vec[i] & trt_grp==trt_grp_vec[j])
    
    # Water year precipitation (Sept-Aug) 
    model_water_year_temp <- lm(abs_cover_pchange ~ water_year_ppt, data=cover_temp)
    summary_water_year_temp <- summary(model_water_year_temp)
    anova_water_year_temp <- Anova(model_water_year_temp, type=3)
    
    out_water_year_temp <- data.frame(trt_grp=trt_grp_vec[j],
                                      PPT_metric="water_yr",
                                      Species = sp_vec[i],
                                      Intercept = summary_water_year_temp$coefficients[1],
                                      Slope = summary_water_year_temp$coefficients[2],
                                      Slope_se = summary_water_year_temp$coefficients[4],
                                      Slope_95lb = confint(model_water_year_temp, 'water_year_ppt')[1],
                                      Slope_95ub = confint(model_water_year_temp, 'water_year_ppt')[2],
                                      df_num=anova_water_year_temp$Df[1],
                                      df_den=anova_water_year_temp$Df[3],
                                      F_val=anova_water_year_temp$F[2],
                                      P_val=anova_water_year_temp$Pr[2],
                                      R2=summary_water_year_temp$r.squared,
                                      adj_R2=summary_water_year_temp$adj.r.squared)
    
    Sens_out_water_year <- rbind(Sens_out_water_year, out_water_year_temp)
    rm(model_water_year_temp, summary_water_year_temp, anova_water_year_temp, out_water_year_temp)
    
    
    rm(cover_temp)
  }
}


### Calculate mean slopes across trt_grp for Table 1.
## Average senstivities across trt_grp
mean_sens_rel <- Sens_out_water_year %>%
  rename(species=Species) %>%
  group_by(species) %>%
  summarise(rel_slope_mean = mean(Slope, na.rm=T),
            rel_slope_se_mean = mean(Slope_se, na.rm=T)) %>%
  mutate(rownum=1:length(levels(factor(Sens_out_water_year$Species))))




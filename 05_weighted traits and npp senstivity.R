### Calculating precipitation senstivity of NPP and compare with community weighted traits - Uses GZTX and LTNPP datasets
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last modified: Nov. 28, 2017

library(tidyverse)
library(ggplot2)
library(ggthemes)

#setwd("T:\\1-People\\Dana Lab\\TRAITS\\PUBLICATIONS\\precip_sensitivity\\")
#setwd("C:\\Users\\wilco\\Desktop\\from ARS comp\\precip_sensitivity\\")
# setwd("C:\\Users\\kwilcox4\\Dropbox\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # work computer
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\")

source("scripts\\02_compile ppt sens_traits_cluster.R")
#rm(list=setdiff(ls(),"full_df"))

source("scripts//01_2_gztx data cleaning.R")
source("scripts//01_3_ltnpp data cleaning.R")
#sort(unique(gztx_cover_raw$Species))


# GZTX --------------------------------------------------------------------

### Identify species to remove (shrubs)
#shrub_vec <- c("ARFR", "ATCA", "CELA", "CHNA", "EREF", "GUSA")

### Calculate species cover averaged over time for each pasture/treatment
gztx_rel_cover_pasture_means <- gztx_cover_raw %>%
#  filter(!Species %in% shrub_vec) %>%
  spread(key=Species, value=abs_cover) %>%
  gather(key=Species, value=abs_cover, -(ID:Quadrat)) %>%
  replace_na(list(abs_cover=0)) %>% 
  group_by(Pasture, Treatment, Species) %>%
  summarise(mean_abs_cover=mean(abs_cover,na.rm=T)) %>%
  mutate(rel_cover=mean_abs_cover/sum(mean_abs_cover))

### Weight traits and senstiivity by species abundances, sum for each pasture
gztx_wtd_ppt_sens <- full_df %>%
  dplyr::select(Species, 
                ppt_slope, 
                LogLDMC, 
                LogLfN, 
                SqrtLfP,
                SqrtHeight,
                SqrtStemSpecDens,
                SqrtPubescence,
                LogSLA,
                LogIndLfArea,
                LogLfThickness,
                LogLfOsmPot) %>%
  right_join(gztx_rel_cover_pasture_means, by="Species") %>%
  mutate(wtd_sens=rel_cover*ppt_slope,
         wtd_LogLDMC=rel_cover*LogLDMC,
         wtd_LogLfN=rel_cover*LogLfN,
         wtd_SqrtLfP=rel_cover*SqrtLfP,
         wtd_SqrtHeight=rel_cover*SqrtHeight,
         wtd_SqrtStemSpecDens=rel_cover*SqrtStemSpecDens,
         wtd_SqrtPubescence=rel_cover*SqrtPubescence,
         wtd_LogSLA=rel_cover*LogSLA,
         wtd_LogIndLfArea=rel_cover*LogIndLfArea,
         wtd_LogLfThickness=rel_cover*LogLfThickness,
         wtd_LogLfOsmPot=rel_cover*LogLfOsmPot) %>%
  group_by(Pasture, Treatment) %>%
  summarise(pasture_ppt_sens = sum( wtd_sens, na.rm=T),
            pasture_LogLDMC = sum( wtd_LogLDMC, na.rm=T),
            pasture_LogLfN = sum( wtd_LogLfN, na.rm=T),
            pasture_SqrtLfP = sum( wtd_SqrtLfP, na.rm=T),
            pasture_SqrtHeight = sum( wtd_SqrtHeight, na.rm=T),
            pasture_SqrtStemSpecDens = sum( wtd_SqrtStemSpecDens, na.rm=T),
            pasture_SqrtPubescence = sum( wtd_SqrtPubescence, na.rm=T),
            pasture_LogSLA = sum( wtd_LogSLA, na.rm=T),
            pasture_LogIndLfArea = sum( wtd_LogIndLfArea, na.rm=T),
            pasture_LogLfThickness = sum( wtd_LogLfThickness, na.rm=T),
            pasture_LogLfOsmPot = sum( wtd_LogLfOsmPot, na.rm=T)  )

gztx_wtd_ppt_sens_long <- gztx_wtd_ppt_sens %>%
  gather(key=trait_name, value=trait_value, -(Pasture:pasture_ppt_sens))


# LTNPP -------------------------------------------------------------------

### Calculate water year precipitation for each year (83-2017)
ppt_water_year <- read.csv("data\\ppt_data\\CPE_dailyprecip_1939-2017Sept.csv") %>%
  mutate(ppt_mm = INCHES*25.4) %>%
  dplyr::select(-INCHES) %>%
  rename(year=YEAR,month=MONTH,day=DAY) %>%
  mutate(water_year = ifelse(month<=8, year, year+1)) %>%
  group_by(water_year) %>%
  summarise(water_year_ppt = sum(ppt_mm,na.rm=T)) %>%
  filter(water_year %in% c(1983:2017)) %>%
  rename(year=water_year)

### Calculate relative abundance for each species in each pasture, averaged across years
ltnpp_rel_cover_pasture_means <- ltnpp_sp_npp %>%
  group_by(site_name, species) %>%
  summarise(mean_anpp=mean(anpp,na.rm=T)) %>%
  mutate(rel_cover=mean_anpp/sum(mean_anpp)) 

### Check the relative cover of the species we have sensitivity values for -- want >90%
sp_vec <- na.omit(full_df)$Species

ltnpp_pcov_of_sp_with_traits <- ltnpp_rel_cover_pasture_means %>%
  mutate(trait_info = ifelse(species %in% sp_vec, "has_sens","no_sens")) %>%
  group_by(site_name, trait_info) %>%
  summarise(sum_relcov=sum(rel_cover))
### All pastures have >95% coverage   
  
### Diagnostics plots (requires only running part of ltnpp_wtd_trt and gztx_std_ppt_sens)
### Looking at relative contribution of species to species-wtd ppt sensitivity values
# ltnpp_wtd_trt$wtd_sp_sens <- with(ltnpp_wtd_trt, ppt_slope*rel_cover)
# ltnpp_wtd_trt <- na.omit(ltnpp_wtd_trt)
# ggplot(ltnpp_wtd_trt, aes(reorder(Species, -wtd_sp_sens), wtd_sp_sens)) +
#   geom_bar(stat="identity") +
#   theme(axis.text.x = element_text(angle=60,hjust=1)) +
#   facet_wrap(~site_name, ncol=1)
### As above, but with gztx
# gztx_wtd_ppt_sens$wtd_sp_sens <- with(gztx_wtd_ppt_sens, ppt_slope*rel_cover)
# ggplot(na.omit(subset(gztx_wtd_ppt_sens, Treatment=="UNUN")), aes(reorder(Species, -wtd_sp_sens), wtd_sp_sens)) +
#   geom_bar(stat="identity") +
#   theme(axis.text.x = element_text(angle=60,hjust=1)) +
#   facet_wrap(~Pasture, ncol=1)
####
# ggplot(subset(ltnpp_pcov_of_sp_with_traits, site_name=="MIDSLOPE"), 
#        aes(reorder(species, -rel_cover), rel_cover, fill=trait_info)) +
#   geom_bar(stat="identity") +
#   theme(axis.text.x=element_text(angle=60, hjust=1))
###

ltnpp_wtd_trt <- full_df %>%
  mutate(Treatment="UNUN") %>%
  dplyr::select(Species,
                Treatment,
                ppt_slope, 
                LogLDMC, 
                LogLfN, 
                SqrtLfP,
                SqrtHeight,
                SqrtStemSpecDens,
                SqrtPubescence,
                LogSLA,
                LogIndLfArea,
                LogLfThickness,
                LogLfOsmPot) %>%
  right_join(ltnpp_rel_cover_pasture_means, by=c("Species"="species")) %>%
  replace_na(list(Treatment="UNUN")) %>% 
  mutate(wtd_sens=rel_cover*ppt_slope,
         wtd_LogLDMC=rel_cover*LogLDMC,
         wtd_LogLfN=rel_cover*LogLfN,
         wtd_SqrtLfP=rel_cover*SqrtLfP,
         wtd_SqrtHeight=rel_cover*SqrtHeight,
         wtd_SqrtStemSpecDens=rel_cover*SqrtStemSpecDens,
         wtd_SqrtPubescence=rel_cover*SqrtPubescence,
         wtd_LogSLA=rel_cover*LogSLA,
         wtd_LogIndLfArea=rel_cover*LogIndLfArea,
         wtd_LogLfThickness=rel_cover*LogLfThickness,
         wtd_LogLfOsmPot=rel_cover*LogLfOsmPot) %>%
  group_by(site_name, Treatment) %>%
  summarise(pasture_ppt_sens=sum(wtd_sens,na.rm=T),
            pasture_LogLDMC=sum(wtd_LogLDMC,na.rm=T),
            pasture_LogLfN=sum(wtd_LogLfN,na.rm=T),
            pasture_SqrtLfP=sum(wtd_SqrtLfP,na.rm=T),
            pasture_SqrtHeight = sum( wtd_SqrtHeight, na.rm=T),
            pasture_SqrtStemSpecDens = sum( wtd_SqrtStemSpecDens, na.rm=T),
            pasture_SqrtPubescence = sum( wtd_SqrtPubescence, na.rm=T),
            pasture_LogSLA=sum(wtd_LogSLA,na.rm=T),
            pasture_LogIndLfArea=sum(wtd_LogIndLfArea,na.rm=T),
            pasture_LogLfThickness=sum(wtd_LogLfThickness,na.rm=T),
            pasture_LogLfOsmPot=sum(wtd_LogLfOsmPot,na.rm=T)) %>%
  rename(Pasture=site_name)

wtd_trt_all <- rbind(gztx_wtd_ppt_sens, ltnpp_wtd_trt)

# Combine GZTX and LTNPP with precipitation data --------------------------------------------------

### calculate precipitaiton senstivity for each pasture and treatment
## precip data ##
ppt_water_year <- read.csv("data\\ppt_data\\CPE_dailyprecip_1939-2017Sept.csv") %>%
  mutate(ppt_mm = INCHES*25.4) %>%
  dplyr::select(-INCHES) %>%
  rename(year=YEAR,month=MONTH,day=DAY) %>%
  mutate(water_year = ifelse(month<=8, year, year+1)) %>%
  group_by(water_year) %>%
  summarise(water_year_ppt = sum(ppt_mm,na.rm=T)) %>%
  filter(water_year %in% c(1983:2017)) %>%
  rename(Year=water_year)

### GZTX ###

### Calculate ANPP (herbaceous) for each pasture in each year (averaged across quadrats)
gztx_npp_data <- read.csv("data\\npp_data\\GZTX_NPP_1992-2015.csv") %>%
  dplyr::select(-NPP, 
                -SpatialLinkToCovDenData, 
                -Weight, 
                -Attribute_LKU, 
                -Derived_Attribute_LKU, 
                -Notes, 
                -FieldMethod_LKU, 
                -LabMethod_LKU,
                -Link_ID) %>%
  rename(Year=YearSampled, Species=species, anpp=WeightPerM2) %>%
#  filter(!Species %in% shrub_vec) %>%
  group_by(Pasture,Treatment, Year, Quadrat) %>%
  summarise(anpp=sum(anpp,na.rm=T)) %>%
  group_by(Pasture, Treatment, Year) %>%
  summarise(anpp=mean(anpp,na.rm=T))

### Combine pasture ANPP data with water-year precipitation 
gztx_npp_ppt <- gztx_npp_data %>%
  left_join(ppt_water_year, by="Year") %>%
  filter(Treatment %in% c("UNUN","GZGZ"))

### LTNPP ### 
ltnpp_npp_ppt <- ltnpp_pasture_npp %>%
  rename(Year=year) %>%
  left_join(ppt_water_year, by="Year") %>%
  mutate(Treatment="UNUN") %>%
  rename(Pasture=site_name) %>%
  dplyr::select(Pasture, Treatment, Year, anpp, water_year_ppt)

### Combine LTNPP with GZTX
npp_ppt <- rbind(gztx_npp_ppt, ltnpp_npp_ppt)
  

#  Calculate sensitivity values -------------------------------------------

pasture_vec <- unique(npp_ppt$Pasture)

sens_out_master_unun <- {}

for(past in 1:length(pasture_vec)){
  
    df_temp <- subset(npp_ppt, Pasture==pasture_vec[past])
    
    ## npp vs Water year precipitation (Sept-Aug)
    # unun -- ungrazed
    unun_model_temp <- lm(anpp ~ water_year_ppt, data=subset(df_temp, Treatment=="UNUN"))
    unun_summary_temp <- summary(unun_model_temp)
    unun_anova_temp <- Anova(unun_model_temp, type=3)
    unun_out_temp <- data.frame(Pasture=pasture_vec[past],
                                PPT_metric="water_yr",
                                Treatment="UNUN",
                                Intercept = unun_model_temp$coefficients[1],
                                Slope = unun_model_temp$coefficients[2],
                                Slope_se = unun_summary_temp$coefficients[4],
                                Slope_95lb = confint(unun_model_temp, 'water_year_ppt')[1],
                                Slope_95ub = confint(unun_model_temp, 'water_year_ppt')[2],
                                df_num=unun_anova_temp$Df[1],
                                df_den=unun_anova_temp$Df[3],
                                F_val=unun_anova_temp$F[2],
                                P_val=unun_anova_temp$Pr[2],
                                R2=unun_summary_temp$r.squared,
                                adj_R2=unun_summary_temp$adj.r.squared)
    sens_out_master_unun <- rbind(sens_out_master_unun, unun_out_temp)
  }

### Grazed pastures only
gzgz_pasture_vec <- unique(subset(npp_ppt, Treatment=="GZGZ")$Pasture)

sens_out_master_gzgz <- {}

for(past in 1:length(gzgz_pasture_vec)){
  
  df_temp <- subset(npp_ppt, Pasture==gzgz_pasture_vec[past])
  
  ## npp vs Water year precipitation (Sept-Aug)
  # gzgz -- mod grazed
  gzgz_model_temp <- lm(anpp ~ water_year_ppt, data=subset(df_temp, Treatment=="GZGZ"))
  gzgz_summary_temp <- summary(gzgz_model_temp)
  gzgz_anova_temp <- Anova(gzgz_model_temp, type=3)
  gzgz_out_temp <- data.frame(Pasture=gzgz_pasture_vec[past],
                              PPT_metric="water_yr",
                              Treatment="GZGZ",
                              Intercept = gzgz_model_temp$coefficients[1],
                              Slope = gzgz_model_temp$coefficients[2],
                              Slope_se = gzgz_summary_temp$coefficients[4],
                              Slope_95lb = confint(gzgz_model_temp, 'water_year_ppt')[1],
                              Slope_95ub = confint(gzgz_model_temp, 'water_year_ppt')[2],
                              df_num=gzgz_anova_temp$Df[1],
                              df_den=gzgz_anova_temp$Df[3],
                              F_val=gzgz_anova_temp$F[2],
                              P_val=gzgz_anova_temp$Pr[2],
                              R2=gzgz_summary_temp$r.squared,
                              adj_R2=gzgz_summary_temp$adj.r.squared)
  sens_out_master_gzgz <- rbind(sens_out_master_gzgz, gzgz_out_temp)
}

sens_out_master <- rbind(sens_out_master_unun, sens_out_master_gzgz)

### standardized slopes
mean_npp <- npp_ppt %>%
  group_by(Pasture, Treatment) %>%
  summarise(mean_anpp=mean(anpp,na.rm=T))

Sens_out_rel_slopes_npp <- sens_out_master %>%
  dplyr::select(Pasture, Treatment, Slope, Slope_se) %>%
  left_join(mean_npp, by=c("Pasture","Treatment")) %>%
  mutate(rel_slope = Slope/mean_anpp,
         rel_slope_se = Slope_se/mean_anpp)

wtd_and_npp_sens <- Sens_out_rel_slopes_npp %>%
  dplyr::select(Pasture, Treatment, Slope, Slope_se, rel_slope, rel_slope_se) %>%
  full_join(wtd_trt_all, by=c("Pasture","Treatment"))
  
wtd_and_npp_sens_long <- wtd_and_npp_sens %>%
  gather(key=trait_name, value=trait_value, -(Pasture:rel_slope_se))



gztx_wtd_sp_sens <- full_df %>%
  dplyr::select(Species,
                ppt_slope,
                LogLDMC,
                LogLfN,
                SqrtLfP,
                SqrtHeight,
                SqrtStemSpecDens,
                SqrtPubescence,
                LogSLA,
                LogIndLfArea,
                LogLfThickness,
                LogLfOsmPot) %>%
  right_join(gztx_rel_cover_pasture_means, by="Species") %>%
  mutate(wtd_sens=rel_cover*ppt_slope,
         wtd_LogLDMC=rel_cover*LogLDMC,
         wtd_LogLfN=rel_cover*LogLfN,
         wtd_SqrtLfP=rel_cover*SqrtLfP,
         wtd_SqrtHeight=rel_cover*SqrtHeight,
         wtd_SqrtStemSpecDens=rel_cover*SqrtStemSpecDens,
         wtd_SqrtPubescence=rel_cover*SqrtPubescence,
         wtd_LogSLA=rel_cover*LogSLA,
         wtd_LogIndLfArea=rel_cover*LogIndLfArea,
         wtd_LogLfThickness=rel_cover*LogLfThickness,
         wtd_LogLfOsmPot=rel_cover*LogLfOsmPot) %>%
  filter(Treatment %in% c("UNUN","GZGZ")) %>%
  mutate(pasture_trt = paste(Pasture, Treatment, sep="_")) %>%
  drop_na()

wtdSLA_plot <- ggplot(gztx_wtd_sp_sens, aes(x=reorder(Species, -wtd_LogSLA ), y=wtd_LogSLA)) +
  geom_bar(stat="identity")+
  facet_wrap(~pasture_trt, scales="fixed")+
  theme_bw() +
  theme(axis.text.x = element_text(angle=45, hjust=1))

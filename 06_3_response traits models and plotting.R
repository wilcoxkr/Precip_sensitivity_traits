### Examines community weighted traits in wet versus dry years - Uses LTGI dataset
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last modified: Feb 16, 2018

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(MASS)
library(lsmeans)

#setwd("T:\\1-People\\Dana Lab\\TRAITS\\PUBLICATIONS\\precip_sensitivity\\")
# setwd("C:\\Users\\kwilcox4\\Dropbox\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # work computer
#setwd("C:\\Users\\wilco\\Dropbox\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # HP laptop
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # HP laptop

source("scripts\\02_compile ppt sens_traits_cluster.R")
#rm(list=setdiff(ls(),"full_df"))

source("scripts//01_1_ltgi data cleaning.R")
#source("scripts//01_2_gztx data cleaning.R")
#source("scripts//01_3_ltnpp data cleaning.R")
#sort(unique(gztx_cover_raw$Species))


emm_options(cov.keep = character(0)) ### allows lsmeans to run as previous package (instead of having predictors treated as 2 level factors)


# ltgi --------------------------------------------------------------------

### Calculate species cover averaged over time for each pasture/treatment
ltgi_rel_cover_annual_means <- ltgi_cover_full %>%
  gather(key=Species, value=abs_cover, -(Year:Plot), -Treatment, -trt_grp) %>%
  replace_na(list(abs_cover=0)) %>% 
  group_by(Pasture, Treatment, Species, Year) %>%
  summarise(mean_abs_cover=mean(abs_cover,na.rm=T)) %>%
  group_by(Pasture, Treatment, Year) %>%
  mutate(rel_cover=mean_abs_cover/sum(mean_abs_cover)) %>%
  replace_na(list(rel_cover=0))

# ggplot(subset(ltgi_rel_cover_annual_means, Pasture=="15E"&Year==2006), aes(x=Species, y=rel_cover)) +
#   geom_bar(stat="identity") +

#   theme(axis.text.x=element_text(angle=45, hjust=1))
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

### Weight traits and senstiivity by species abundances, sum for each pasture
ltgi_wtd_traits <- full_df %>%
  dplyr::select(Species, 
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
  right_join(ltgi_rel_cover_annual_means, by="Species") %>%
  mutate(wtd_LogLDMC=rel_cover*LogLDMC,
         wtd_LogLfN=rel_cover*LogLfN,
         wtd_SqrtLfP=rel_cover*SqrtLfP,
         wtd_SqrtHeight=rel_cover*SqrtHeight,
         wtd_SqrtStemSpecDens=rel_cover*SqrtStemSpecDens,
         wtd_SqrtPubescence=rel_cover*SqrtPubescence,
         wtd_LogSLA=rel_cover*LogSLA,
         wtd_LogIndLfArea=rel_cover*LogIndLfArea,
         wtd_LogLfThickness=rel_cover*LogLfThickness,
         wtd_LogLfOsmPot=rel_cover*LogLfOsmPot) %>%
  group_by(Pasture, Treatment, Year) %>%
  summarise(pasture_LogLDMC = sum( wtd_LogLDMC, na.rm=T),
            pasture_LogLfN = sum( wtd_LogLfN, na.rm=T),
            pasture_SqrtLfP = sum( wtd_SqrtLfP, na.rm=T),
            pasture_SqrtHeight = sum( wtd_SqrtHeight, na.rm=T),
            pasture_SqrtStemSpecDens = sum( wtd_SqrtStemSpecDens, na.rm=T),
            pasture_SqrtPubescence = sum( wtd_SqrtPubescence, na.rm=T),
            pasture_LogSLA = sum( wtd_LogSLA, na.rm=T),
            pasture_LogIndLfArea = sum( wtd_LogIndLfArea, na.rm=T),
            pasture_LogLfThickness = sum( wtd_LogLfThickness, na.rm=T),
            pasture_LogLfOsmPot = sum( wtd_LogLfOsmPot, na.rm=T)  )

ltgi_wtd_traits_long <- ltgi_wtd_traits %>%
  mutate(trt_grp = ifelse(Treatment %in% c("UN","LG"), "ULG", "MHG")) %>% # lumps ungrazed with light grazed & mod with heavy grazed
  gather(key=trait_name, value=trait_value, -(Pasture:Year), -trt_grp) %>%
  group_by(trt_grp, Year, trait_name) %>%
  summarize(trait_value=mean(trait_value))


### Combine with water year precipitatoin
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

  ### merge with weighted trait annual data
ltgi_wtd_traits_ppt_long <- ltgi_wtd_traits_long %>%
  rename(year=Year) %>%
  left_join(ppt_water_year, by="year")

### Scale data for each trait
### Put data in long form and normalize trait values
ltgi_scaled_long <- ltgi_wtd_traits_ppt_long %>%
  group_by(trait_name) %>%
  mutate(
    water_year_ppt_norm = scale(water_year_ppt),
    trait_value_norm = scale(trait_value)
  )

ltgi_scaled_long$trt_grp <- factor(ltgi_scaled_long$trt_grp, levels=c("ULG","MHG"))
ggplot(ltgi_scaled_long, aes(water_year_ppt, trait_value_norm, col=trt_grp)) +
  geom_point()+
  geom_smooth(method="lm",se=F)+
  facet_wrap(~trait_name, scales="free") +
  theme_few()

### Run model loop to generate slopes and se's for each life form and functional group
trait_vec <- unique(ltgi_scaled_long$trait_name) ## vector of traits to cycle through

ppt_summary_out <- {} ## set up data frame to fill with model output
ppt_anova_out <- {}
ppt_lstrends_out <- {}

for(trait in 1:length(trait_vec)){ # loop through traits
  
  ltgi_df_temp <- ltgi_scaled_long %>% # create temporary data frame for one trait at a time
    filter(trait_name == trait_vec[trait])
  
  # run models
  ppt_model_temp <- lm(trait_value_norm ~ water_year_ppt_norm*trt_grp, data=ltgi_df_temp)
  
  ppt_model_summary <- summary(ppt_model_temp)
  ppt_model_anova <- Anova(ppt_model_temp, type=3)
  ppt_lstrends <- summary(emtrends(ppt_model_temp, "trt_grp", var="water_year_ppt_norm"))

  # write model output
  ppt_summary_out_temp <- data.frame( Regression = "PPT-Tcomm",
                                      Out_type = "summary",
                              Trait_name = trait_vec[trait],
                              Source = rownames(ppt_model_summary$coefficients),
                              Slope = ppt_model_summary$coefficients[,1],
                              SE = ppt_model_summary$coefficients[,2],
                              CI_lower = confint(ppt_model_temp)[,1],
                              CI_upper = confint(ppt_model_temp)[,2],
                              t_value = ppt_model_summary$coefficients[,3],
                              P_value = ppt_model_summary$coefficients[,4])
  
  ppt_anova_out_temp <- data.frame( Regression = "PPT-Tcomm",
                                  Out_type = "Anova",
                                  Trait_name = trait_vec[trait],
                                  Source = rownames(ppt_model_anova),
                                  Df=ppt_model_anova$Df,
                                  F_value=ppt_model_anova$`F value`,
                                  P_value=ppt_model_anova$`Pr(>F)`)

  ppt_lstrends_out_temp <- data.frame(
                                      Regression = "PPT-Tcomm",
                                      Out_type = "lstrends",
                                      Trait_name=trait_vec[trait],
                                     trt_grp=ppt_lstrends$trt_grp,
                                     Df=ppt_lstrends$df,
                                     Slope=ppt_lstrends$water_year_ppt_norm.trend,
                                     SE=ppt_lstrends$SE,
                                     CI_lower = ppt_lstrends$lower.CL,
                                     CI_upper = ppt_lstrends$upper.CL
                                     )
  
  ppt_summary_out <- rbind(ppt_summary_out, ppt_summary_out_temp)
  ppt_anova_out <- rbind(ppt_anova_out, ppt_anova_out_temp)
  ppt_lstrends_out <- rbind(ppt_lstrends_out, ppt_lstrends_out_temp)
  }


### PLOTTTING CODE ###
ggplot(ppt_lstrends_out, aes(x=Trait_name, y=Slope, ymin=CI_lower, ymax=CI_upper, col=trt_grp)) +
  geom_point(size=3, position=position_dodge(width=0.2)) +
  geom_errorbar(position=position_dodge(width=0.2), width=0) +
  geom_hline(yintercept=0, lty=2) +
  theme_classic() +
  theme(axis.text.x = element_text(angle=45, hjust=1))


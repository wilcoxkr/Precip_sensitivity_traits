### Plotting species-weighted versus npp sensitivity
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### last updated: Dec. 12, 2017


#setwd("T:\\1-People\\Dana Lab\\TRAITS\\PUBLICATIONS\\precip_sensitivity\\")
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\")
library(tidyverse)
library(ggthemes)
library(gridExtra)

source("scripts//05_weighted traits and npp senstivity.R")
rm(list=setdiff(ls(),c("wtd_and_npp_sens","npp_ppt")))


### Main figure


## regression line is weighted but points are not sized accordingly - WITH GRAZED PLOTS
npp_sp_plot <- ggplot(subset(wtd_and_npp_sens,
                               Pasture!= "19"), 
                      aes(x=pasture_ppt_sens, y=Slope)) +
  geom_point(size=3) +
  geom_smooth(method="lm",se=F, col="black",
              mapping = aes(weight = 1/Slope_se), 
              show.legend=F, 
              data=subset(wtd_and_npp_sens, !Pasture %in% c("19","ESA","ESA_no_ARFR"))) +
  theme_few() +
  xlab("Weighted species-level sensitivity") +
  ylab("ANPP sensitivity (g m-2 mm-1)")

pdf("figures//npp vs wtd sp sens_with grazed_inclESAnoARFR_Feb2017.pdf",height=3, width=3, useDingbats = F)
print(npp_sp_plot)
dev.off()

### Main model(s)

### Create subsets for models and weighting
df_ungrazed_no19 <- wtd_and_npp_sens %>%
  dplyr::select(Pasture, Treatment, Slope, Slope_se, pasture_ppt_sens) %>%
  filter(Treatment=="UNUN") %>%
  filter(!Pasture %in% c(19))

df_ungrazed_no19noESA <- wtd_and_npp_sens %>%
  dplyr::select(Pasture, Treatment, Slope, Slope_se, pasture_ppt_sens) %>%
  filter(Treatment=="UNUN") %>%
  filter(!Pasture %in% c(19, "ESA"))

df_no19 <- wtd_and_npp_sens %>%
  dplyr::select(Pasture, Treatment, Slope, Slope_se, pasture_ppt_sens) %>%
  filter(!Pasture %in% c(19))

df_no19noESA <- wtd_and_npp_sens %>%
  dplyr::select(Pasture, Treatment, rel_slope, rel_slope_se, Slope, Slope_se, pasture_ppt_sens) %>%
  filter(!Pasture %in% c(19, "ESA"))

## ungrazed no pasture 19 model
weights_ungrazed_no19 <- 1/df_ungrazed_no19$Slope_se
lm_ungrazed_no19 <- lm(Slope ~ pasture_ppt_sens,
              weights=weights_ungrazed_no19,
              data=df_ungrazed_no19)
summary(lm_ungrazed_no19)
## ungrazed no pasture 19 or ESA in model
weights_ungrazed_no19noESA <- 1/df_ungrazed_no19noESA$Slope_se
lm_ungrazed_no19noESA <- lm(Slope ~ pasture_ppt_sens,
              weights=weights_ungrazed_no19noESA,
              data=df_ungrazed_no19noESA)
summary(lm_ungrazed_no19noESA)

## ungrazed and grazed no pasture 19 model
weights_no19 <- 1/df_no19$Slope_se
lm_no19 <- lm(Slope ~ pasture_ppt_sens,
                       weights=weights_no19,
                       data=df_no19)
summary(lm_no19)

## ungrazed no pasture 19 or ESA in model
weights_no19noESA <- 1/df_no19noESA$Slope_se
lm_no19noESA <- lm(Slope ~ pasture_ppt_sens,
                            weights=weights_no19noESA,
                            data=df_no19noESA)
summary(lm_no19noESA)


### relative anpp slope
## ungrazed no pasture 19 or ESA in model
weights_no19noESA <- 1/df_no19noESA$rel_slope_se
lm_no19noESA <- lm(rel_slope ~ pasture_ppt_sens,
                   weights=weights_no19noESA,
                   data=df_no19noESA)
summary(lm_no19noESA)

npp_sp_plot_rel <- ggplot(subset(wtd_and_npp_sens,
                             Pasture!= "19"), 
                      aes(x=pasture_ppt_sens, y=rel_slope, label=Pasture)) +
#  geom_point(size=3) +
  geom_text() +
  geom_smooth(method="lm",se=F, col="black",
              mapping = aes(weight = 1/rel_slope_se), 
              show.legend=F, 
              data=subset(wtd_and_npp_sens, !Pasture %in% c("19","ESA","ESA_no_ARFR"))) +
  theme_few() +
  xlab("Weighted species-level sensitivity") +
  ylab("ANPP sensitivity (g m-2 mm-1)")

pdf("figures//npp vs wtd sp sens_with grazed_rel anpp slope_Feb2017.pdf",height=3, width=3, useDingbats = F)
print(npp_sp_plot_rel)
dev.off()

png("figures//npp vs wtd sp sens_with grazed_rel anpp slope_Feb2017.png",height=3, width=3, units="in", res=600)
print(npp_sp_plot_rel)
dev.off()



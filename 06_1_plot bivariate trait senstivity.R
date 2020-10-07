### Plotting bivariate trait-sensitivity relationships
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last modified: Nov. 30, 2017


### Set up workspace
rm(list=ls())
#setwd("T:\\1-People\\Dana Lab\\TRAITS\\PUBLICATIONS\\precip_sensitivity\\")
#setwd("C:\\Users\\wilco\\Desktop\\from ARS comp\\precip_sensitivity\\")
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # HP laptop

library(tidyverse)
library(ggthemes)

### Read in data
source("scripts//02_compile ppt sens_traits_cluster.R")
rm(list=setdiff(ls(),"full_df"))

### Define focal traits
trait_vec_all <- c(
  "LogSLA",
  "LogLfN",
  "SqrtLfP",
  "SqrtStemSpecDens",
  "LogLfOsmPot",
  "LogLfThickness",
  "SqrtPubescence",
  "LogLDMC",
  "LogIndLfArea",
  "SqrtHeight"
)


### Put data in long form with only focal traits
full_df_long <- full_df %>%
  drop_na() %>%
  gather(key=trait_name,value=trait_value,-(Species:life_form)) %>%
  filter(trait_name %in% trait_vec_all) %>%
  mutate(trait_name=factor(trait_name, levels=c(trait_vec_all) ))

### plot differentiated by life form
full_plot <- ggplot(full_df_long, aes(x=trait_value, y=ppt_slope, size=1/ppt_slope_se)) +
  geom_point(aes(x=trait_value, y=ppt_slope, col=life_form), alpha=0.7) +
  geom_smooth(method = "lm", se=F, mapping = aes(weight = 1/ppt_slope_se), show.legend=F)  + 
  theme_few() + 
  xlab("Trait value") +
  ylab("Precipitation senstivity (cover change (%)/(ppt change (mm) x mean cover (%)))") +
  facet_wrap(~trait_name, scales="free")

pdf("figures//trait senstivity biplots_by life form.pdf", width=10.5, height=6,  useDingbats = F)
print(full_plot)
dev.off()

### Plot d13C
d13c_df <- full_df %>%
  drop_na() %>%
  gather(key=trait_name,value=trait_value,-(Species:life_form)) %>%
  filter(trait_name=="LogLf13C") %>%
  mutate(photo_path= ifelse(trait_value<1, "C4", "C3")) %>%
  mutate(photo_path=factor(photo_path, levels=c("C4","C3")))

# shapes
d13c_plot_1 <- 
  ggplot(d13c_df, aes(x=trait_value, y=ppt_slope, size=1/ppt_slope_se, show.legend=F)) +
  geom_point(aes(x=trait_value, y=ppt_slope, col=life_form)) +
  geom_smooth(method = "lm", se=F, mapping = aes(weight = 1/ppt_slope_se), show.legend=F)  + 
  theme_few() + 
  theme(legend.position="none") +
  xlab("Trait value") +
  ylab("Sens") +
  facet_wrap(~photo_path, scales="free_x")

pdf("figures//trait senstivity biplots_d13C_color.pdf", width=2.449, height=2.191,  useDingbats = F)
  print(d13c_plot_1)
dev.off()
  


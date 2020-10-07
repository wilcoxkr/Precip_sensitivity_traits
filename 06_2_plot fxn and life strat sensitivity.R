### Plotting bivariate trait-sensitivity relationships
### Author: Kevin Wilcox (wilcoxkr@gmail.com)
###
### Last modified: Nov. 30, 2017


### Set up workspace
rm(list=ls())
#setwd("T:\\1-People\\Dana Lab\\TRAITS\\PUBLICATIONS\\precip_sensitivity\\")
setwd("C:\\Users\\wilco\\OneDrive - University of Wyoming\\Cross_workstation_workspace\\Current projects\\precip_sensitivity\\") # HP laptop

library(tidyverse)
library(ggplot2)
library(ggthemes)
library(lsmeans)

### Read in data
source("scripts//02_compile ppt sens_traits_cluster.R")
rm(list=setdiff(ls(),"full_df"))


### Calculate mean and standard error for each functional group and life strategy
se=function(x){sd(x[!is.na(x)])/sqrt(length(x[!is.na(x)]))}

fxn_grp_stats <- full_df %>%
  drop_na() %>%
  group_by(fxn_grp) %>%
  summarise_at(vars(ppt_slope),funs(mean, se)) %>%
  rename(group=fxn_grp)

life_form_stats <- full_df %>%
  drop_na() %>%
  group_by(life_form) %>%
  summarise_at(vars(ppt_slope),funs(mean, se)) %>%
  rename(group=life_form)

all_stats <- rbind(fxn_grp_stats, life_form_stats)

# group_plot <- ggplot(all_stats, aes(x=group, y=mean, ymin=mean-se, ymax=mean+se)) +
#   geom_hline(yintercept=0, col="dark grey") +
#   geom_bar(stat="identity", col="black", fill="lightgrey") +
#   geom_errorbar(width=.1) +
#   theme_few() +
#   ylab("Precipitation senstivity (dCover/(dPPT*muCover))") +
#   xlab("Functional group or Life strategy")
# 
# 
# pdf("figures//fxn life sensitivity barplot_Nov2017.pdf", width=4, height=3,  useDingbats = F)
# print(group_plot)
# dev.off()

### Run model and calculate LSMEANS

fxn_grp_weights <- 1/full_df$ppt_slope_se
fxn_grp_model <- lm(ppt_slope ~ fxn_grp + life_form,
                    weights=fxn_grp_weights,
                    data=full_df)
Anova(fxn_grp_model,type=3)
lsmeans_out <- summary(lsmeans(fxn_grp_model, list("fxn_grp", "life_form"), contr = "pairwise"))

fxn_lsmeans <- data.frame(group=lsmeans_out$`lsmeans of fxn_grp`$fxn_grp,
                          lsmean=lsmeans_out$`lsmeans of fxn_grp`$lsmean,
                          se=lsmeans_out$`lsmeans of fxn_grp`$SE,
                          ci_lower=lsmeans_out$`lsmeans of fxn_grp`$lower.CL,
                          ci_upper=lsmeans_out$`lsmeans of fxn_grp`$upper.CL)
life_lsmeans <- data.frame(group=lsmeans_out$`lsmeans of life_form`$life_form,
                          lsmean=lsmeans_out$`lsmeans of life_form`$lsmean,
                          se=lsmeans_out$`lsmeans of life_form`$SE,
                          ci_lower=lsmeans_out$`lsmeans of life_form`$lower.CL,
                          ci_upper=lsmeans_out$`lsmeans of life_form`$upper.CL)
all_lsmeans <- rbind(fxn_lsmeans, life_lsmeans) %>%
  mutate(group=factor(group, levels=c("Forb", "Sshrub", "Graminoid", "Annual", "Perennial")))


lsmeans_plot <- ggplot(all_lsmeans, aes(x=group, y=lsmean, ymin=lsmean-se, ymax=lsmean+se)) +
  geom_hline(yintercept=0, col="dark grey") +
  geom_col(col="black", fill="lightgrey") +
  geom_errorbar(width=.1) +
  theme_few() +
  ylab("Precipitation senstivity (dCover/(dPPT*muCover))") +
  xlab("Functional group or Life strategy")

pdf("figures//fxn life sensitivity lsmeans barplot_Nov2017.pdf", width=4, height=3,  useDingbats = F)
print(lsmeans_plot)
dev.off()




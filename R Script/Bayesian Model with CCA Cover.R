### Bayesian Model with CCA Cover ####

library(scatterplot3d); library(reshape); library (reshape2); library (data.table); 
library(RColorBrewer); library(dplyr); library (tidyr); library (plyr); library (ggplot2); 
require (vegan); require (goeveg);   require  (stats); require (dynRB); require (matrixStats); 
require(plyr); require (stringr);require (reshape2); require(gridExtra)
require (betapart); library (car); library (MASS); library (glmm); library (glmnet); 
library (glmmTMB); library (lme4); library(randomcoloR)
require (ggpubr); library (cowplot); library (patchwork); library (scales); 
library(viridisLite) ; library(patchwork); library(fishualize);library(tidyverse) 
require (DHARMa); require (MuMIn); require (merTools)

# ***Mid-domain effect (there's no mid-domain effect) ??? ####
require(devtools)
library(reshape2)

# ***Geospatial data ??? ####
require (geodist)

# rm (list = ls())
# getwd()
# setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R")

### Set working directory ####
BC.cover <- read.csv("BC.cover.csv", header=T, sep=",")
# if want to contain only "Crustose Coralline Algae"
CCA.cover <- subset (BC.cover, MajorCategory == "Crustose Coralline Algae")
CCA.cover <- aggregate (Cover ~ Region + Site + Depth + Transect + Quadrat + MajorCategory + Genus, CCA.cover , sum)

CCA.cover_df <- CCA.cover

# Modify Coral Cover from per picture 100% to per site 100%
CCA.cover <- ddply(CCA.cover, ~ Region + Site + Depth + MajorCategory + Genus, function(x){c(Cover = sum(x$Cover))})
CCA.cover <- CCA.cover %>%
  mutate(Pictures = case_when(
    Depth == 5 ~ 105,  # 105 pictures in 5 m in each site
    Depth == 15 ~ 105, # 105 pictures in 15 m in each site
    Depth == 30 ~ 63)) # 63 pictures in 30 m in each site
CCA.cover <- CCA.cover %>%
  group_by(Region, Site, Depth, MajorCategory, Genus) %>%
  reframe(Cover = Cover/Pictures)

CCA_data <- CCA.cover


### Coral cover profile ####
# Plot considering the effect of Region and Site to see the effect of depth
# Calculate sum for coral cover in function of site and depth for a quick interpretation and plot
CCA_cover <- aggregate (Cover ~ Region + Site + Depth, CCA_data , sum)
# For keeping quadrats - perhaps unnecessary
CCA.cover_df <- aggregate (Cover ~ Region + Site + Quadrat +  Depth, CCA.cover_df , sum)

# Transform depth as a Qualitative variable  
CCA_cover$Depth <- as.factor(as.character(CCA_cover$Depth))
CCA_cover$Depth = factor(CCA_cover$Depth,levels = c ("5", "15", "30"))

# ggplot with Sites (Location)
ggplot(CCA_cover, aes(x=Depth, y=Cover)) + 
  geom_boxplot() + geom_point(aes (),size = 1) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=4) + 
  theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

## no-mid-domain effect.


### Mixed-Model for showing cover decreases with depth in all sites ####
# Just follow the part of modify the data
CCA_cover2 <- CCA_cover
CCA_cover2$Depth <- as.numeric (as.character(CCA_cover2$Depth))
CCA_cover2$Site <- as.character(CCA_cover2$Site)
CCA_cover2$Location <- paste(CCA_cover2$Site, "_", CCA_cover2$Depth)

### Use a beta distribution to avoid problem of normality ####
# Just follow the part of modify the data
CCA_cover2$Cover_Beta <- CCA_cover2$Cover/100

CCA_cover2$Tot_Points <- 25
CCA_cover2$Coral_points <- (CCA_cover2$Cover * CCA_cover2$Tot_Points) / 100
CCA_cover2$Coral_points  <- round(CCA_cover2$Coral_points,0) #calculate how many point fall on coral to the nearest integer
CCA_cover2$NonCoral_points <- CCA_cover2$Tot_Points - CCA_cover2$Coral_points 

CCA_cover2$Proportion = CCA_cover2$Coral_points / (CCA_cover2$Coral_points + CCA_cover2$NonCoral_points)
CCA_cover2$Proportion [CCA_cover2$Proportion == 0] <-  .001 # Otherwise beta distribution with 0 it does not work


# The best solution is to use Bayesian modelling. ####
library(brms);library('rstan'); library("stam");library("parallel"); library ("performance")
# ***I think no prior here!!! ####
Binomial_CCA.Cover_model <- brm(Coral_points | trials(Tot_Points) ~  Depth + (1 | Site),
                               data = CCA_cover2, family = binomial(),  # prior = my_priors,
                               control = list(adapt_delta = 0.9, max_treedepth = 11),
                               iter = 4000, warmup = 1000, chains = 2, cores = 2) 
save(Binomial_CCA.Cover_model, file="/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Binomial_CCA.Cover_model_Diversity.RData")
load("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Binomial_CCA.Cover_model_Diversity.RData")
# they test: Binomial_Cover_model, Beta_Cover_model, Beta_Cover_model_zero, and Student_model

# The binomial is better than the beta model ####
summary (Binomial_CCA.Cover_model)
bayes_R2(Binomial_CCA.Cover_model)

# The best is to use the Binomial_Cover_Model ####
plot(Binomial_CCA.Cover_model)
# "pp_check" help to evaluate how well the Bayesian model fits the data by comparing the observed data with data simulated from the posterior distribution of the model parameters
pp_check(Binomial_CCA.Cover_model, type = "scatter_avg") # Not structured data
bayes_R2(Binomial_CCA.Cover_model) # 
r2_bayes(Binomial_CCA.Cover_model)

CCA.me_null <- conditional_effects(Binomial_CCA.Cover_model, nsamples = 1000, probs = c(0.05, 0.95), spaghetti = F) # Default is 0.95
plot(CCA.me_null, ask = FALSE, points = F) # Probability scale!

# ***don't understand the part of posterior predict! ####
# Work with posterior predict ####
# Create the reference dataframe - newdata (Here for all depths)
Depth <- unique (Binomial_CCA.Cover_model$data$Depth)
Site <- unique (Binomial_CCA.Cover_model$data$Site)
Coral_points <- unique (Binomial_CCA.Cover_model$data$Coral_points)
Tot_Points <- unique (Binomial_CCA.Cover_model$data$Tot_Points)
#
CCA.ref_data <- crossing(Depth, Site, Tot_Points)

# Forest plots with the posteriors
CCA.fitted_values <- posterior_epred(Binomial_CCA.Cover_model, newdata = CCA.ref_data, re_formula = 'Coral_points | trials(Tot_Points) ~  Depth + (1 | Site)')

# Number of rows equals to CCA.ref_data and number of dimensions is equal to (Site*Depths=27)
# str (CCA.fitted_values)
# dim (CCA.fitted_values)

# Necessary to traspose
CCA.fitted_values <- t(CCA.fitted_values)

# Create combination of CCA.ref_data
CCA.ref_data_fitted <- cbind (CCA.ref_data [c(1,2)],CCA.fitted_values)

# Multiple columns of predictions into a single one. 
CCA.ref_data_fitted <- melt (CCA.ref_data_fitted, id.vars = c ("Depth", "Site"), na.rm = F, measure.vars = c(3:6002), value.name = c("Posterior_Prob"))
CCA.ref_data_fitted$Depth = factor (CCA.ref_data_fitted$Depth, levels = c ("30", "15","5"))


# Summary depth and island among all iterations 

CCA.summary <- ddply(CCA.ref_data_fitted, .(Depth,Site), summarize, Post_Mean=mean(Posterior_Prob), Post_Sd = sd(Posterior_Prob), Post_se=sd(Posterior_Prob) / sqrt(length(Posterior_Prob)), 
                    Post_Margin.error = qt(p=0.05/2, df=length (Posterior_Prob)-1,lower.tail=F) * Post_se)

CCA_cover2 <- merge (CCA_cover2, CCA.summary)



# Extract the number of points from cover / Same as done before
CCA_cover2$Tot_Points <- 25
CCA_cover2$Coral_points <- (CCA_cover2$Cover * CCA_cover2$Tot_Points) / 100
# Round the points
CCA_cover2$Coral_points <-round(CCA_cover2$Coral_points,0)
CCA_cover2$NonCoral_points <-abs (CCA_cover2$Coral_points - CCA_cover2$Tot_Points)


# Transform the posterior predict of binomial from points proportions (out of 25) to coral cover (%)
CCA_cover2$Post_Mean_Cover <- (CCA_cover2$Post_Mean * 100) / CCA_cover2$Tot_Points
CCA_cover2$Post_Sd_Cover <- (CCA_cover2$Post_Sd * 100) / CCA_cover2$Tot_Points
CCA_cover2$Post_se_Cover <- (CCA_cover2$Post_se * 100) / CCA_cover2$Tot_Points
CCA_cover2$Post_Margin.error_Cover <- (CCA_cover2$Post_Margin.error * 100) / CCA_cover2$Tot_Points



#################################### new #############################

# Measure mean and standard error and confidence intervals from raw cover values

CCA2_cover <- CCA_cover2
CCA2.summary <- ddply(CCA2_cover, .(Depth), summarize, Mean=mean(Cover), Cover_se=sd(Cover) / sqrt(length(Cover)), 
                     Margin.error = qt(p=0.05/2, df=length (Cover)-1,lower.tail=F) * Cover_se)

# Confidence Interval manually would be
# alpha = 0.05
# degrees.freedom = 16 - 1
# t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
# print(t.score)
# margin.error <- t.score * summary$Cover_se
# # Confidence Interval

# Combine dataframes
CCA2_cover <- merge (CCA2_cover,CCA2.summary)

# Make a single value for each depth of Post Bayesian
CCA2.summary2 <- ddply(CCA2_cover, .(Depth), summarize, Post_Cover_Depth=mean(Post_Mean_Cover), Post_Cover_sd_Depth = sd (Post_Mean_Cover))
CCA2_cover <- merge (CCA2_cover,CCA2.summary2)


# Make the plot again
CCA2_cover$Depth <- as.numeric (as.character(CCA2_cover$Depth))



# Plot with Bayesian Binomial model (with standard deviation of the model and with coral cover mean values with CI

CCA.Fig_1A <- ggplot(CCA2_cover, aes(x=Depth, y=Cover)) +
  # geom_point(aes(fill = Island, colour = Island),shape = 21, size = 0.8)  + 
  geom_errorbar(aes(ymin = Mean - Margin.error, ymax = Mean + Margin.error), color="black", size=1, width=2) +
  geom_line(aes(y=Post_Cover_Depth), size=1, linetype="dashed", alpha = 1) +
  geom_line(aes(y=Post_Cover_Depth + Post_Cover_sd_Depth), size=0.3, linetype="dotted",alpha = 0.8) +
  geom_line(aes(y=Post_Cover_Depth - Post_Cover_sd_Depth), size=0.3, linetype="dotted",alpha = 0.8) +
  scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
  geom_point(aes(y=Mean), shape=21, fill="white", size=4) +  
  scale_x_continuous(name ="Depth (m)", limits=c(3,35), breaks = c(5,15,30)) +
  scale_y_continuous(name ="CCA cover (%)", limits=c(-5,55), breaks = c(0,10,20,30,40,50)) +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "none") 
CCA.Fig_1A


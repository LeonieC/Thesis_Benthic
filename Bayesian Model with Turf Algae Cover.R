### Bayesian Model with Turf Algae Cover ####

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
# if want to contain only "Turf Algae" we defined
TA.cover <- subset (BC.cover, MajorCategory == "Turf Algae")
TA.cover <- aggregate (Cover ~ Region + Site + Depth + Transect + Quadrat + MajorCategory + Genus, TA.cover , sum)

TA.cover_df <- TA.cover

# Modify Coral Cover from per picture 100% to per site 100%
TA.cover <- ddply(TA.cover, ~ Region + Site + Depth + MajorCategory + Genus, function(x){c(Cover = sum(x$Cover))})
TA.cover <- TA.cover %>%
  mutate(Pictures = case_when(
    Depth == 5 ~ 105,  # 105 pictures in 5 m in each site
    Depth == 15 ~ 105, # 105 pictures in 15 m in each site
    Depth == 30 ~ 63)) # 63 pictures in 30 m in each site
TA.cover <- TA.cover %>%
  group_by(Region, Site, Depth, MajorCategory, Genus) %>%
  reframe(Cover = Cover/Pictures)

TA_data <- TA.cover


### Coral cover profile ####
# Plot considering the effect of Region and Site to see the effect of depth
# Calculate sum for coral cover in function of site and depth for a quick interpretation and plot
TA_cover <- aggregate (Cover ~ Region + Site + Depth, TA_data , sum)
# For keeping quadrats - perhaps unnecessary
TA.cover_df <- aggregate (Cover ~ Region + Site + Quadrat +  Depth, TA.cover_df , sum)

# Transform depth as a Qualitative variable  
TA_cover$Depth <- as.factor(as.character(TA_cover$Depth))
TA_cover$Depth = factor(TA_cover$Depth,levels = c ("5", "15", "30"))

# ggplot with Sites (Location)
ggplot(TA_cover, aes(x=Depth, y=Cover)) + 
  geom_boxplot() + geom_point(aes (),size = 1) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=4) + 
  theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

## no-mid-domain effect.


### Mixed-Model for showing cover decreases with depth in all sites ####
# Just follow the part of modify the data
TA_cover2 <- TA_cover
TA_cover2$Depth <- as.numeric (as.character(TA_cover2$Depth))
TA_cover2$Site <- as.character(TA_cover2$Site)
TA_cover2$Location <- paste(TA_cover2$Site, "_", TA_cover2$Depth)

### Use a beta distribution to avoid problem of normality ####
# Just follow the part of modify the data
TA_cover2$Cover_Beta <- TA_cover2$Cover/100

TA_cover2$Tot_Points <- 25
TA_cover2$Coral_points <- (TA_cover2$Cover * TA_cover2$Tot_Points) / 100
TA_cover2$Coral_points  <- round(TA_cover2$Coral_points,0) #calculate how many point fall on coral to the nearest integer
TA_cover2$NonCoral_points <- TA_cover2$Tot_Points - TA_cover2$Coral_points 

TA_cover2$Proportion = TA_cover2$Coral_points / (TA_cover2$Coral_points + TA_cover2$NonCoral_points)
TA_cover2$Proportion [TA_cover2$Proportion == 0] <-  .001 # Otherwise beta distribution with 0 it does not work


# The best solution is to use Bayesian modelling. ####
library(brms);library('rstan'); library("stam");library("parallel"); library ("performance")
# ***I think no prior here!!! ####
Binomial_TA.Cover_model <- brm(Coral_points | trials(Tot_Points) ~  Depth + (1 | Site),
                               data = TA_cover2, family = binomial(),  # prior = my_priors,
                               control = list(adapt_delta = 0.9, max_treedepth = 11),
                               iter = 4000, warmup = 1000, chains = 2, cores = 2) 
save(Binomial_TA.Cover_model, file="/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Binomial_TA.Cover_model_Diversity.RData")
load("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Binomial_TA.Cover_model_Diversity.RData")
# they test: Binomial_Cover_model, Beta_Cover_model, Beta_Cover_model_zero, and Student_model

# The binomial is better than the beta model ####
summary (Binomial_TA.Cover_model)
bayes_R2(Binomial_TA.Cover_model)

# The best is to use the Binomial_Cover_Model ####
plot(Binomial_TA.Cover_model)
# "pp_check" help to evaluate how well the Bayesian model fits the data by comparing the observed data with data simulated from the posterior distribution of the model parameters
pp_check(Binomial_TA.Cover_model, type = "scatter_avg") # Not structured data
bayes_R2(Binomial_TA.Cover_model) # 
r2_bayes(Binomial_TA.Cover_model)

TA.me_null <- conditional_effects(Binomial_TA.Cover_model, nsamples = 1000, probs = c(0.05, 0.95), spaghetti = F) # Default is 0.95
plot(TA.me_null, ask = FALSE, points = F) # Probability scale!

# ***don't understand the part of posterior predict! ####
# Work with posterior predict ####
# Create the reference dataframe - newdata (Here for all depths)
Depth <- unique (Binomial_TA.Cover_model$data$Depth)
Site <- unique (Binomial_TA.Cover_model$data$Site)
Coral_points <- unique (Binomial_TA.Cover_model$data$Coral_points)
Tot_Points <- unique (Binomial_TA.Cover_model$data$Tot_Points)
#
TA.ref_data <- crossing(Depth, Site, Tot_Points)

# Forest plots with the posteriors
TA.fitted_values <- posterior_epred(Binomial_TA.Cover_model, newdata = TA.ref_data, re_formula = 'Coral_points | trials(Tot_Points) ~  Depth + (1 | Site)')

# Number of rows equals to TA.ref_data and number of dimensions is equal to (Site*Depths=27)
# str (TA.fitted_values)
# dim (TA.fitted_values)

# Necessary to traspose
TA.fitted_values <- t(TA.fitted_values)

# Create combination of TA.ref_data
TA.ref_data_fitted <- cbind (TA.ref_data [c(1,2)],TA.fitted_values)

# Multiple columns of predictions into a single one. 
TA.ref_data_fitted <- melt (TA.ref_data_fitted, id.vars = c ("Depth", "Site"), na.rm = F, measure.vars = c(3:6002), value.name = c("Posterior_Prob"))
TA.ref_data_fitted$Depth = factor (TA.ref_data_fitted$Depth, levels = c ("30", "15","5"))


# Summary depth and island among all iterations 

TA.summary <- ddply(TA.ref_data_fitted, .(Depth,Site), summarize, Post_Mean=mean(Posterior_Prob), Post_Sd = sd(Posterior_Prob), Post_se=sd(Posterior_Prob) / sqrt(length(Posterior_Prob)), 
                    Post_Margin.error = qt(p=0.05/2, df=length (Posterior_Prob)-1,lower.tail=F) * Post_se)

TA_cover2 <- merge (TA_cover2, TA.summary)



# Extract the number of points from cover / Same as done before
TA_cover2$Tot_Points <- 25
TA_cover2$Coral_points <- (TA_cover2$Cover * TA_cover2$Tot_Points) / 100
# Round the points
TA_cover2$Coral_points <-round(TA_cover2$Coral_points,0)
TA_cover2$NonCoral_points <-abs (TA_cover2$Coral_points - TA_cover2$Tot_Points)


# Transform the posterior predict of binomial from points proportions (out of 25) to coral cover (%)
TA_cover2$Post_Mean_Cover <- (TA_cover2$Post_Mean * 100) / TA_cover2$Tot_Points
TA_cover2$Post_Sd_Cover <- (TA_cover2$Post_Sd * 100) / TA_cover2$Tot_Points
TA_cover2$Post_se_Cover <- (TA_cover2$Post_se * 100) / TA_cover2$Tot_Points
TA_cover2$Post_Margin.error_Cover <- (TA_cover2$Post_Margin.error * 100) / TA_cover2$Tot_Points



#################################### new #############################

# Measure mean and standard error and confidence intervals from raw cover values

TA2_cover <- TA_cover2
TA2.summary <- ddply(TA2_cover, .(Depth), summarize, Mean=mean(Cover), Cover_se=sd(Cover) / sqrt(length(Cover)), 
                     Margin.error = qt(p=0.05/2, df=length (Cover)-1,lower.tail=F) * Cover_se)

# Confidence Interval manually would be
# alpha = 0.05
# degrees.freedom = 16 - 1
# t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
# print(t.score)
# margin.error <- t.score * summary$Cover_se
# # Confidence Interval

# Combine dataframes
TA2_cover <- merge (TA2_cover,TA2.summary)

# Make a single value for each depth of Post Bayesian
TA2.summary2 <- ddply(TA2_cover, .(Depth), summarize, Post_Cover_Depth=mean(Post_Mean_Cover), Post_Cover_sd_Depth = sd (Post_Mean_Cover))
TA2_cover <- merge (TA2_cover,TA2.summary2)


# Make the plot again
TA2_cover$Depth <- as.numeric (as.character(TA2_cover$Depth))



# Plot with Bayesian Binomial model (with standard deviation of the model and with coral cover mean values with CI

TA.Fig_1A <- ggplot(TA2_cover, aes(x=Depth, y=Cover)) +
  # geom_point(aes(fill = Island, colour = Island),shape = 21, size = 0.8)  + 
  geom_errorbar(aes(ymin = Mean - Margin.error, ymax = Mean + Margin.error), color="black", size=1, width=2) +
  geom_line(aes(y=Post_Cover_Depth), size=1, linetype="dashed", alpha = 1) +
  geom_line(aes(y=Post_Cover_Depth + Post_Cover_sd_Depth), size=0.3, linetype="dotted",alpha = 0.8) +
  geom_line(aes(y=Post_Cover_Depth - Post_Cover_sd_Depth), size=0.3, linetype="dotted",alpha = 0.8) +
  scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
  geom_point(aes(y=Mean), shape=21, fill="white", size=4) +  
  scale_x_continuous(name ="Depth (m)", limits=c(3,35), breaks = c(5,15,30)) +
  scale_y_continuous(name ="Turf Algae cover (%)", limits=c(15,75), breaks = c(20,30,40,50,60,70)) +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "none") 
TA.Fig_1A


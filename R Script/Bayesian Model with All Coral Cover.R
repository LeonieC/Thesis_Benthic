### Bayesian Model with All Coral Cover ####

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
# if want to contain all "Coral" we defined
Coral.cover <- subset(BC.cover, MajorCategory %in% c("Scleractinian", "Antipatharian", "Octocoral", "Gorgonian"))
Coral.cover <- aggregate (Cover ~ Region + Site + Depth + Transect + Quadrat + MajorCategory + Genus, Coral.cover , sum)

Coral.cover_df <- Coral.cover

# Modify Coral Cover from per picture 100% to per site 100%
Coral.cover <- ddply(Coral.cover, ~ Region + Site + Depth + MajorCategory + Genus, function(x){c(Cover = sum(x$Cover))})
Coral.cover <- Coral.cover %>%
  mutate(Pictures = case_when(
    Depth == 5 ~ 105,  # 105 pictures in 5 m in each site
    Depth == 15 ~ 105, # 105 pictures in 15 m in each site
    Depth == 30 ~ 63)) # 63 pictures in 30 m in each site
Coral.cover <- Coral.cover %>%
  group_by(Region, Site, Depth, MajorCategory, Genus) %>%
  reframe(Cover = Cover/Pictures)

Coral_data <- Coral.cover


### Coral cover profile ####
# Plot considering the effect of Region and Site to see the effect of depth
# Calculate sum for coral cover in function of site and depth for a quick interpretation and plot
Coral_cover <- aggregate (Cover ~ Region + Site + Depth, Coral_data , sum)
# For keeping quadrats - perhaps unnecessary
Coral.cover_df <- aggregate (Cover ~ Region + Site + Quadrat +  Depth, Coral.cover_df , sum)

# Transform depth as a Qualitative variable  
Coral_cover$Depth <- as.factor(as.character(Coral_cover$Depth))
Coral_cover$Depth = factor(Coral_cover$Depth,levels = c ("5", "15", "30"))

# ggplot with Sites (Location)
ggplot(Coral_cover, aes(x=Depth, y=Cover)) + 
  geom_boxplot() + geom_point(aes (),size = 1) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=4) + 
  theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))

## no-mid-domain effect.


### Mixed-Model for showing cover decreases with depth in all sites ####
# Just follow the part of modify the data
Coral_cover2 <- Coral_cover
Coral_cover2$Depth <- as.numeric (as.character(Coral_cover2$Depth))
Coral_cover2$Site <- as.character(Coral_cover2$Site)
Coral_cover2$Location <- paste(Coral_cover2$Site, "_", Coral_cover2$Depth)

### Use a beta distribution to avoid problem of normality ####
# Just follow the part of modify the data
Coral_cover2$Cover_Beta <- Coral_cover2$Cover/100

Coral_cover2$Tot_Points <- 25
Coral_cover2$Coral_points <- (Coral_cover2$Cover * Coral_cover2$Tot_Points) / 100
Coral_cover2$Coral_points  <- round(Coral_cover2$Coral_points,0) #calculate how many point fall on coral to the nearest integer
Coral_cover2$NonCoral_points <- Coral_cover2$Tot_Points - Coral_cover2$Coral_points 

Coral_cover2$Proportion = Coral_cover2$Coral_points / (Coral_cover2$Coral_points + Coral_cover2$NonCoral_points)
Coral_cover2$Proportion [Coral_cover2$Proportion == 0] <-  .001 # Otherwise beta distribution with 0 it does not work


# The best solution is to use Bayesian modelling. ####
library(brms);library('rstan'); library("stam");library("parallel"); library ("performance")
# ***I think no prior here!!! ####
Binomial_Coral.cover_model <- brm(Coral_points | trials(Tot_Points) ~  Depth + (1 | Site),
                               data = Coral_cover2, family = binomial(),  # prior = my_priors,
                               control = list(adapt_delta = 0.9, max_treedepth = 11),
                               iter = 4000, warmup = 1000, chains = 2, cores = 2) 
save(Binomial_Coral.cover_model, file="/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Binomial_Coral.cover_model_Diversity.RData")
load("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Binomial_Coral.cover_model_Diversity.RData")
# they test: Binomial_Cover_model, Beta_Cover_model, Beta_Cover_model_zero, and Student_model

# The binomial is better than the beta model ####
summary (Binomial_Coral.cover_model)
bayes_R2(Binomial_Coral.cover_model)

# The best is to use the Binomial_Cover_Model ####
plot(Binomial_Coral.cover_model)
# "pp_check" help to evaluate how well the Bayesian model fits the data by comparing the observed data with data simulated from the posterior distribution of the model parameters
pp_check(Binomial_Coral.cover_model, type = "scatter_avg") # Not structured data
bayes_R2(Binomial_Coral.cover_model) # 
r2_bayes(Binomial_Coral.cover_model)

Coral.me_null <- conditional_effects(Binomial_Coral.cover_model, nsamples = 1000, probs = c(0.05, 0.95), spaghetti = F) # Default is 0.95
plot(Coral.me_null, ask = FALSE, points = F) # Probability scale!

# ***don't understand the part of posterior predict! ####
# Work with posterior predict ####
# Create the reference dataframe - newdata (Here for all depths)
Depth <- unique (Binomial_Coral.cover_model$data$Depth)
Site <- unique (Binomial_Coral.cover_model$data$Site)
Coral_points <- unique (Binomial_Coral.cover_model$data$Coral_points)
Tot_Points <- unique (Binomial_Coral.cover_model$data$Tot_Points)
#
Coral.ref_data <- crossing(Depth, Site, Tot_Points)

# Forest plots with the posteriors
Coral.fitted_values <- posterior_epred(Binomial_Coral.cover_model, newdata = Coral.ref_data, re_formula = 'Coral_points | trials(Tot_Points) ~  Depth + (1 | Site)')

# Number of rows equals to Coral.ref_data and number of dimensions is equal to (Site*Depths=27)
# str (Coral.fitted_values)
# dim (Coral.fitted_values)

# Necessary to traspose
Coral.fitted_values <- t(Coral.fitted_values)

# Create combination of Coral.ref_data
Coral.ref_data_fitted <- cbind (Coral.ref_data [c(1,2)],Coral.fitted_values)

# Multiple columns of predictions into a single one. 
Coral.ref_data_fitted <- melt (Coral.ref_data_fitted, id.vars = c ("Depth", "Site"), na.rm = F, measure.vars = c(3:6002), value.name = c("Posterior_Prob"))
Coral.ref_data_fitted$Depth = factor (Coral.ref_data_fitted$Depth, levels = c ("30", "15","5"))


# Summary depth and island among all iterations 

Coral.summary <- ddply(Coral.ref_data_fitted, .(Depth,Site), summarize, Post_Mean=mean(Posterior_Prob), Post_Sd = sd(Posterior_Prob), Post_se=sd(Posterior_Prob) / sqrt(length(Posterior_Prob)), 
                    Post_Margin.error = qt(p=0.05/2, df=length (Posterior_Prob)-1,lower.tail=F) * Post_se)

Coral_cover2 <- merge (Coral_cover2, Coral.summary)



# Extract the number of points from cover / Same as done before
Coral_cover2$Tot_Points <- 25
Coral_cover2$Coral_points <- (Coral_cover2$Cover * Coral_cover2$Tot_Points) / 100
# Round the points
Coral_cover2$Coral_points <-round(Coral_cover2$Coral_points,0)
Coral_cover2$NonCoral_points <-abs (Coral_cover2$Coral_points - Coral_cover2$Tot_Points)


# Transform the posterior predict of binomial from points proportions (out of 25) to coral cover (%)
Coral_cover2$Post_Mean_Cover <- (Coral_cover2$Post_Mean * 100) / Coral_cover2$Tot_Points
Coral_cover2$Post_Sd_Cover <- (Coral_cover2$Post_Sd * 100) / Coral_cover2$Tot_Points
Coral_cover2$Post_se_Cover <- (Coral_cover2$Post_se * 100) / Coral_cover2$Tot_Points
Coral_cover2$Post_Margin.error_Cover <- (Coral_cover2$Post_Margin.error * 100) / Coral_cover2$Tot_Points



#################################### new #############################

# Measure mean and standard error and confidence intervals from raw cover values

Coral2_cover <- Coral_cover2
Coral2.summary <- ddply(Coral2_cover, .(Depth), summarize, Mean=mean(Cover), Cover_se=sd(Cover) / sqrt(length(Cover)), 
                     Margin.error = qt(p=0.05/2, df=length (Cover)-1,lower.tail=F) * Cover_se)

# Confidence Interval manually would be
# alpha = 0.05
# degrees.freedom = 16 - 1
# t.score = qt(p=alpha/2, df=degrees.freedom,lower.tail=F)
# print(t.score)
# margin.error <- t.score * summary$Cover_se
# # Confidence Interval

# Combine dataframes
Coral2_cover <- merge (Coral2_cover,Coral2.summary)

# Make a single value for each depth of Post Bayesian
Coral2.summary2 <- ddply(Coral2_cover, .(Depth), summarize, Post_Cover_Depth=mean(Post_Mean_Cover), Post_Cover_sd_Depth = sd (Post_Mean_Cover))
Coral2_cover <- merge (Coral2_cover,Coral2.summary2)


# Make the plot again
Coral2_cover$Depth <- as.numeric (as.character(Coral2_cover$Depth))



# Plot with Bayesian Binomial model (with standard deviation of the model and with coral cover mean values with CI

Coral.Fig_1A <- ggplot(Coral2_cover, aes(x=Depth, y=Cover)) +
  # geom_point(aes(fill = Island, colour = Island),shape = 21, size = 0.8)  + 
  geom_errorbar(aes(ymin = Mean - Margin.error, ymax = Mean + Margin.error), color="black", size=1, width=2) +
  geom_line(aes(y=Post_Cover_Depth), size=1, linetype="dashed", alpha = 1) +
  geom_line(aes(y=Post_Cover_Depth + Post_Cover_sd_Depth), size=0.3, linetype="dotted",alpha = 0.8) +
  geom_line(aes(y=Post_Cover_Depth - Post_Cover_sd_Depth), size=0.3, linetype="dotted",alpha = 0.8) +
  scale_fill_manual(values = cols) + scale_color_manual(values = cols) +
  geom_point(aes(y=Mean), shape=21, fill="white", size=4) +  
  scale_x_continuous(name ="Depth (m)", limits=c(3,35), breaks = c(5,15,30)) +
  scale_y_continuous(name ="Coral cover (%)", limits=c(-5,55), breaks = c(0,10,20,30,40,50)) +
  theme_classic() + theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
                          axis.text = element_text(size=10, colour="black"),
                          axis.title = element_text(size=11, face="bold", colour="black"), legend.position = "none") 
Coral.Fig_1A


### Bayesian Model -- All Coral Cover with Light Attenuation ####

# Data manipulation and aggregation / Visualization
library(plyr); library(dplyr); library(ggplot2)
# Bayesian modeling
library(brms); library(rstan); library(stam); library(parallel)
# Model performance evaluation
library(performance)


rm (list = ls())
getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file")

### Set up the working data ####
BC.cover <- read.csv("BC.cover.csv", header=T, sep=",")
# to contain all "Corals" we defined
Coral.cover <- subset(BC.cover, MajorCategory %in% c("Black Corals", "Gorgonian Corals", "Soft Corals", "Stony Corals"))
Coral.cover <- aggregate (Cover ~ Region + Site + Depth + Transect + Quadrat + MajorCategory + OTUs, Coral.cover , sum)

# Modify Coral Cover from per picture 100% to per location (Site + Depth) 100%
Coral.cover <- ddply(Coral.cover, ~ Region + Site + Depth + MajorCategory + OTUs, function(x){c(Cover = sum(x$Cover))})
Coral.cover <- Coral.cover %>%
  mutate(Pictures = case_when(
    Depth == 5 ~ 105,  # 105 pictures in 5 m in each site
    Depth == 15 ~ 105, # 105 pictures in 15 m in each site
    Depth == 30 ~ 63)) # 63 pictures in 30 m in each site
Coral.cover <- Coral.cover %>%
  group_by(Region, Site, Depth, MajorCategory, OTUs) %>%
  reframe(Cover = Cover/Pictures)

Coral_data <- Coral.cover


### Coral cover profile ####
# Plot considering the effect of Region and Site to see the effect of depth
# Calculate sum for coral cover in function of site and depth for a quick interpretation and plot
Coral_cover <- aggregate (Cover ~ Region + Site + Depth, Coral_data , sum)

# Transform depth as a Qualitative variable  
Coral_cover$Depth <- as.factor(as.character(Coral_cover$Depth))
Coral_cover$Depth = factor(Coral_cover$Depth,levels = c ("5", "15", "30"))

# ggplot with Locations
ggplot(Coral_cover, aes(x=Depth, y=Cover)) + 
  geom_boxplot() + geom_point(aes (),size = 1) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=4) + 
  theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))
## mid-domain effect???



### Modify the data for Bayesian model with PAR ####
Coral_cover2 <- Coral_cover
Coral_cover2$Depth <- as.numeric (as.character(Coral_cover2$Depth))
Coral_cover2$Site <- as.character(Coral_cover2$Site)
Coral_cover2$Location <- paste(Coral_cover2$Site, "_", Coral_cover2$Depth)
Coral_cover2$RegDep <- paste(Coral_cover2$Region, "_", Coral_cover2$Depth)
# Add PAR value (use PAR to replace Depth)
PAR.Location <- read.csv("PAR.Location.csv", header=T, sep=",")
Coral_cover2 <- Coral_cover2 %>%
  left_join(dplyr::select(PAR.Location, Region, Site, Depth, PAR.D.Mean, PAR.S.Mean, Light.Attenuation), by = c("Region", "Site", "Depth"))

# Use a beta distribution to avoid problem of normality
Coral_cover2$Cover_Beta <- Coral_cover2$Cover/100

Coral_cover2$Tot_Points <- 25
Coral_cover2$Coral_points <- (Coral_cover2$Cover * Coral_cover2$Tot_Points) / 100
Coral_cover2$Coral_points  <- round(Coral_cover2$Coral_points,0) #calculate how many point fall on coral to the nearest integer
Coral_cover2$NonCoral_points <- Coral_cover2$Tot_Points - Coral_cover2$Coral_points 

Coral_cover2$Proportion = Coral_cover2$Coral_points / (Coral_cover2$Coral_points + Coral_cover2$NonCoral_points)
Coral_cover2$Proportion [Coral_cover2$Proportion == 0] <-  .001 # Otherwise beta distribution with 0 it does not work



# Bayesian modelling ####
Binomial_Corals.CoverLightA_model <- brm(Coral_points | trials(Tot_Points) ~  Light.Attenuation + (1 | Site),
                               data = Coral_cover2, family = binomial(),  # prior = my_priors,
                               control = list(adapt_delta = 0.9, max_treedepth = 11),
                               iter = 4000, warmup = 1000, chains = 2, cores = 2) 
# save(Binomial_Corals.CoverLightA_model, file="/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Binomial_Corals.CoverLightA_model_Diversity.RData")
load("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file/Binomial_Corals.CoverLightA_model_Diversity.RData")


# Summary & Performance evaluation ####
summary (Binomial_Corals.CoverLightA_model)
bayes_R2(Binomial_Corals.CoverLightA_model)
plot(Binomial_Corals.CoverLightA_model)

# "pp_check" help to evaluate how well the Bayesian model fits the data by comparing the observed data with data simulated from the posterior distribution of the model parameters
pp_check(Binomial_Corals.CoverLightA_model, type = "scatter_avg") # Not structured data
bayes_R2(Binomial_Corals.CoverLightA_model) # 
r2_bayes(Binomial_Corals.CoverLightA_model)

Coral.me_null <- conditional_effects(Binomial_Corals.CoverLightA_model, nsamples = 1000, probs = c(0.05, 0.95), spaghetti = F) # Default is 0.95
plot(Coral.me_null, ask = FALSE, points = F) # Probability scale!


# Plot Bayesian Model ####
# Extract the data from Coral.me_null
coral_effects_df <- as.data.frame(Coral.me_null$Light.Attenuation)
# Predictive plot with confidence interval
bayes_r2_val <- bayes_R2(Binomial_Corals.CoverLightA_model)[1]

# Point Color by Region
Coral_cover2$Region <- factor(Coral_cover2$Region, levels = c("North", "Green Island", "Xiaoliuqiu"))

plot1 <- ggplot(coral_effects_df, aes(x = Light.Attenuation, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "lightgray", alpha = 0.4) +  # Confidence interval
  geom_line(linewidth = 1) +  # Main line
  geom_point(data = Coral_cover2, aes(x = Light.Attenuation, y = Cover_Beta, color = Region), size = 2) + # Real data points
  theme_minimal() +
  labs(x = "% of Surface PAR", y = "Cover (%)", color = "Region") +  # Add color legend label
  scale_y_continuous(limits = c(0, 1),
                     breaks = seq(0, 1, 0.2),
                     labels = seq(0, 100, 20)) +
  scale_x_continuous(limits = c(0, 0.36),
                     breaks = seq(0, 0.36, 0.1),
                     labels = seq(0, 36, 10)) +
  scale_color_manual(
    values = c("#98c1d9", "#ffc857", "#adc178"),
    labels = c("NT", "GI", "XLQ")) +
  theme(
    axis.title = element_text(size = 14, face = "bold"),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(color = "#ced4da"),
    panel.grid.minor = element_line(color = "#e9ecef"),
    axis.title.x = element_blank(),      # Removes x-axis title
    axis.text.x = element_blank(),        # Removes x-axis text/ticks
    legend.position = "none"  # Removes the legend
  ) +
  annotate("text", x = 0.3, y = 0.85, label = paste("Bayesian R² : ", round(bayes_r2_val, 3)), size = 3) +
  annotate("text", x = 0.3, y = 0.95, label = paste("All Corals"), size = 4, fontface = "bold")

plot1



### Bayesian Model -- CCA Cover with Light Attenuation ####

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
# to contain only "Crustose Coralline Algae" we defined
CCA.cover <- subset(BC.cover, MajorCategory %in% c("Crustose Coralline Algae"))
CCA.cover <- aggregate (Cover ~ Region + Site + Depth + Transect + Quadrat + MajorCategory + OTUs, CCA.cover , sum)

# Modify Coral Cover from per picture 100% to per location (Site + Depth) 100%
CCA.cover <- ddply(CCA.cover, ~ Region + Site + Depth + MajorCategory + OTUs, function(x){c(Cover = sum(x$Cover))})
CCA.cover <- CCA.cover %>%
  mutate(Pictures = case_when(
    Depth == 5 ~ 105,  # 105 pictures in 5 m in each site
    Depth == 15 ~ 105, # 105 pictures in 15 m in each site
    Depth == 30 ~ 63)) # 63 pictures in 30 m in each site
CCA.cover <- CCA.cover %>%
  group_by(Region, Site, Depth, MajorCategory, OTUs) %>%
  reframe(Cover = Cover/Pictures)

CCA_data <- CCA.cover


### Coral cover profile ####
# Plot considering the effect of Region and Site to see the effect of depth
# Calculate sum for coral cover in function of site and depth for a quick interpretation and plot
CCA_cover <- aggregate (Cover ~ Region + Site + Depth, CCA_data , sum)

# Transform depth as a Qualitative variable  
CCA_cover$Depth <- as.factor(as.character(CCA_cover$Depth))
CCA_cover$Depth = factor(CCA_cover$Depth,levels = c ("5", "15", "30"))

# ggplot with Locations
ggplot(CCA_cover, aes(x=Depth, y=Cover)) + 
  geom_boxplot() + geom_point(aes (),size = 1) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=4) + 
  theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))
## mid-domain effect???



### Modify the data for Bayesian model with PAR ####
CCA_cover2 <- CCA_cover
CCA_cover2$Depth <- as.numeric (as.character(CCA_cover2$Depth))
CCA_cover2$Site <- as.character(CCA_cover2$Site)
CCA_cover2$Location <- paste(CCA_cover2$Site, "_", CCA_cover2$Depth)
CCA_cover2$RegDep <- paste(CCA_cover2$Region, "_", CCA_cover2$Depth)
# Add PAR value (use PAR to replace Depth)
PAR.Location <- read.csv("PAR.Location.csv", header=T, sep=",")
CCA_cover2 <- CCA_cover2 %>%
  left_join(dplyr::select(PAR.Location, Region, Site, Depth, PAR.D.Mean, PAR.S.Mean, Light.Attenuation), by = c("Region", "Site", "Depth"))

# Use a beta distribution to avoid problem of normality
CCA_cover2$Cover_Beta <- CCA_cover2$Cover/100

CCA_cover2$Tot_Points <- 25
CCA_cover2$Coral_points <- (CCA_cover2$Cover * CCA_cover2$Tot_Points) / 100
CCA_cover2$Coral_points  <- round(CCA_cover2$Coral_points,0) #calculate how many point fall on coral to the nearest integer
CCA_cover2$NonCoral_points <- CCA_cover2$Tot_Points - CCA_cover2$Coral_points 

CCA_cover2$Proportion = CCA_cover2$Coral_points / (CCA_cover2$Coral_points + CCA_cover2$NonCoral_points)
CCA_cover2$Proportion [CCA_cover2$Proportion == 0] <-  .001 # Otherwise beta distribution with 0 it does not work



# Bayesian modelling ####
Binomial_CCA.CoverLightA_model <- brm(Coral_points | trials(Tot_Points) ~  Light.Attenuation + (1 | Site),
                               data = CCA_cover2, family = binomial(),  # prior = my_priors,
                               control = list(adapt_delta = 0.9, max_treedepth = 11),
                               iter = 4000, warmup = 1000, chains = 2, cores = 2) 
# save(Binomial_CCA.CoverLightA_model, file="/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Binomial_CCA.CoverLightA_model_Diversity.RData")
load("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file/Binomial_CCA.CoverLightA_model_Diversity.RData")


# Summary & Performance evaluation ####
summary (Binomial_CCA.CoverLightA_model)
bayes_R2(Binomial_CCA.CoverLightA_model)
plot(Binomial_CCA.CoverLightA_model)

# "pp_check" help to evaluate how well the Bayesian model fits the data by comparing the observed data with data simulated from the posterior distribution of the model parameters
pp_check(Binomial_CCA.CoverLightA_model, type = "scatter_avg") # Not structured data
bayes_R2(Binomial_CCA.CoverLightA_model) # 
r2_bayes(Binomial_CCA.CoverLightA_model)

CCA.me_null <- conditional_effects(Binomial_CCA.CoverLightA_model, nsamples = 1000, probs = c(0.05, 0.95), spaghetti = F) # Default is 0.95
plot(CCA.me_null, ask = FALSE, points = F) # Probability scale!


# Plot Bayesian Model ####
# Extract the data from CCA.me_null
CCA_effects_df <- as.data.frame(CCA.me_null$Light.Attenuation)
# Predictive plot with confidence interval
bayes_r2_val <- bayes_R2(Binomial_CCA.CoverLightA_model)[1]

# Point Color by Region
CCA_cover2$Region <- factor(CCA_cover2$Region, levels = c("North", "Green Island", "Xiaoliuqiu"))

plot6 <- ggplot(CCA_effects_df, aes(x = Light.Attenuation, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "lightgray", alpha = 0.4) +  # Confidence interval
  geom_line(linewidth = 1) +  # Main line
  geom_point(data = CCA_cover2, aes(x = Light.Attenuation, y = Cover_Beta, color = Region), size = 2) +  # Real data points
  theme_minimal() +
  labs(x = "% of Surface PAR", y = "Cover (%)") +
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
    axis.title.y = element_blank(),      # Removes y-axis title
    axis.text.y = element_blank(),       # Removes y-axis text/ticks
    legend.title = element_text(size = 11.5, face = "bold"),  # Adjust legend title
    legend.text = element_text(size = 9.5)                 # Adjust legend text
  ) +
  annotate("text", x = 0.3, y = 0.85, label = paste("Bayesian R² : ", round(bayes_r2_val, 3)), size = 3) +
  annotate("text", x = 0.27, y = 0.95, label = paste("Crustose Coralline Algae"), size = 4, fontface = "bold")

plot6


# * for putting all Bayesian plots together ####
library(patchwork)
# plotA <- plot1 + plot4
# plotB <- plot2 + plot5
# plotC <- plot3 + plot6
# plotA / plotB / plotC

# Combine all plots and add legend at the bottom
# combined_plot <- (plotA / plotB / plotC) + plot_layout(guides = "collect")
# combined_plot



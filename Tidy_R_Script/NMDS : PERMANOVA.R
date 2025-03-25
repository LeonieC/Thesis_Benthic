### NMDS / PERMANOVA of BC & CA ####
## beta quantitative => NMDS + bray curtis + permanova + simper (indicative species)
library(dplyr); library(tidyr); library(tibble)
library(vegan); library(labdsv); library(RVAideMemoire)
# Load the pairwiseAdonis package if not already done
# install.packages("devtools")
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

rm (list = ls())
getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file")

### Benthic Community NMDS---in OTMU level ####
BC.only <- read.csv("BC.only.csv", header=T, sep=",")

OTMU_bc.sum <- BC.only %>% 
  group_by(Region, Location, Depth, MajorCategory, OTMUs) %>% 
  summarise(Number = length(OTMUs))

OTMU_bc.sum$Region <- factor(OTMU_bc.sum$Region, 
                            levels = c('North', 'Green Island','Xiaoliuqiu'))

simi.bc <- OTMU_bc.sum %>%
  mutate(Region = factor(Region, levels = unique(Region)), 
         Location = factor(Location, levels = unique(Location))) %>%
  ungroup() %>%
  group_by(Location,OTMUs) %>% summarize(Abundance = sum(Number))%>%
  pivot_wider(names_from = c(OTMUs), values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  column_to_rownames(var = "Location")

rowSums(simi.bc)
simi.bc <- as.matrix(simi.bc)
OTMU_bc.cover <- proportions(simi.bc, 1)
head(OTMU_bc.cover)
rowSums(OTMU_bc.cover)
OTMU_bc.cover <- as.data.frame(OTMU_bc.cover)

#create Factor table
Factor <- OTMU_bc.sum
Location<-row.names(OTMU_bc.cover)
head(Location)
RegionDepth<-c("North15","North30","North5","North15","North30","North5","North15","North30","North5", 
               "Green Island15","Green Island30","Green Island5","Green Island15","Green Island30","Green Island5","Green Island15","Green Island30","Green Island5", 
               "Xiaoliuqiu15","Xiaoliuqiu30","Xiaoliuqiu5","Xiaoliuqiu15","Xiaoliuqiu30","Xiaoliuqiu5","Xiaoliuqiu15","Xiaoliuqiu30","Xiaoliuqiu5")
Region<-c("North","North","North","North","North","North","North","North","North", 
          "Green Island","Green Island","Green Island","Green Island","Green Island","Green Island","Green Island","Green Island","Green Island", 
          "Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu")
Depth<-c("15","30","5","15","30","5","15","30","5",
         "15","30","5","15","30","5","15","30","5",
         "15","30","5","15","30","5","15","30","5")

Factor<-as.data.frame(cbind(Location,RegionDepth,Region, Depth))
head(Factor)
Factor$RegionDepth<-as.factor(Factor$RegionDepth)
Factor$Location<-as.factor(Factor$Location)
Factor$Region<-as.factor(Factor$Region)
Factor$Depth<-as.factor(Factor$Depth)
#transform cover data due to many 0
OTMU_bc.cover_hell<-decostand(OTMU_bc.cover,"hellinger")

#Meta NMDS
nmds1 <- metaMDS(OTMU_bc.cover, "bray", type='n')

## *indicator species analysis of BC ####
siv <- indval(OTMU_bc.cover_hell, Factor$Region)
gr <- siv$maxcls[siv$pval <= 0.05]
iv <- siv$indcls[siv$pval <= 0.05]
pv <- siv$pval[siv$pval <= 0.05]
fidg <- data.frame(group = gr, indval = iv, pvalue = pv)
fidg <- fidg[order(fidg$group, -fidg$indval), ]
# write.csv(fidg, 'indicator species of BC.csv',row.names = T)

# Add group names to fidg
group_levels <- levels(Factor$Region) # make sure the group fit to right Region (check levels)
fidg$group_name <- factor(fidg$group, labels = group_levels)

# Extract species scores for NMDS plot
species_scores <- as.data.frame(vegan::scores(nmds1, display = "species"))
species_scores$species <- rownames(species_scores)

# Adjust plotting area and margins
xlim_range <- c(-1, 2)  # Set custom x-axis limits
ylim_range <- range(c(species_scores$NMDS2, species_scores$NMDS2))  # Adjust y-axis limits as needed
par(mar = c(4, 4, 1, 1))  # Increase margins to provide space for text

# Plot NMDS of BC ####
plot(nmds1, display = "sites", type = 'n', xlab = "NMDS1", ylab = "NMDS2", xlim = xlim_range, ylim = ylim_range)
# Add hulls for different regions
ordihull(nmds1, groups = Factor$Region, col = c("#ffc857", "#98c1d9", "#adc178"), lwd = 2)
# Add points for different regions and depths
points(nmds1, "sites", col = "#C2F1F2", pch = 17, select = Factor$RegionDepth=="North5")
points(nmds1, "sites", col = "#98c1d9", pch = 15, select = Factor$RegionDepth=="North15")
points(nmds1, "sites", col = "#274c77", pch = 16, select = Factor$RegionDepth=="North30")
points(nmds1, "sites", col = "#FFEE8C", pch = 17, select = Factor$RegionDepth=="Green Island5")
points(nmds1, "sites", col = "#ffc857", pch = 15, select = Factor$RegionDepth=="Green Island15")
points(nmds1, "sites", col = "#ff9f1c", pch = 16, select = Factor$RegionDepth=="Green Island30")
points(nmds1, "sites", col = "#EAF291", pch = 17, select = Factor$RegionDepth=="Xiaoliuqiu5")
points(nmds1, "sites", col = "#adc178", pch = 15, select = Factor$RegionDepth=="Xiaoliuqiu15")
points(nmds1, "sites", col = "#709775", pch = 16, select = Factor$RegionDepth=="Xiaoliuqiu30")

## *stress of BC NMDS => nmds1$stress ####
# Add stress value
mtext(paste("Stress:", round(nmds1$stress, 2)), side = 1, line = 3, cex = 0.8, adj = 1)

# Add significant indicator species to NMDS plot
for (i in 1:nrow(fidg)) {
  species <- rownames(fidg)[i]
  if (species %in% rownames(species_scores)) {
    text(species_scores[species, 1], species_scores[species, 2], labels = species, cex = 0.5, pos = 4, col = "#979dac")
  }
}
# ***can't add legend :( => i want two legend: Depth symbol + Region color ####


## PERMANOVA of BC ####
# PERMANOVA -> test difference between Region? between Depth? or between Region*Depth?

#Region
PermRegion<-adonis2(OTMU_bc.cover~Region, data=Factor, method="bray", permutations=999)
pairwise.perm.manova(vegdist(OTMU_bc.cover,"bray"), Factor$Region, nperm=9999)

#Depth
PermDepth<-adonis2(OTMU_bc.cover~Depth, data=Factor, method="bray", permutations=999)
pairwise.perm.manova(vegdist(OTMU_bc.cover,"bray"), Factor$Depth, nperm=9999)

#Depth*Region (interaction between D & R)
Factor$DepthRegion <- interaction(Factor$Depth, Factor$Region)
PermDepthRegion<-adonis2(OTMU_bc.cover~Depth*Region, data=Factor, method="bray", permutations=999)
pairwise.perm.manova(vegdist(OTMU_bc.cover,"bray"), Factor$DepthRegion, nperm=9999)


BC.Perm.R <- PermRegion
BC.Perm.D <- PermDepth
BC.Perm.DR <- PermDepthRegion
# write.csv(BC.Perm.R, 'BC.Perm.R.csv',row.names = T)
# write.csv(BC.Perm.D, 'BC.Perm.D.csv',row.names = T)
# write.csv(BC.Perm.DR, 'BC.Perm.DR.csv',row.names = T)

## pairwise PERMANOVA of BC ####
BC.pairPermR <- pairwise.adonis(vegdist(OTMU_bc.cover, "bray"), factors = Factor$Region, perm = 9999)
BC.pairPermD <- pairwise.adonis(vegdist(OTMU_bc.cover, "bray"), factors = Factor$Depth, perm = 9999)
BC.pairPermDR <- pairwise.adonis(vegdist(OTMU_bc.cover, "bray"), factors = Factor$DepthRegion, perm = 9999)
# write.csv(BC.pairPermR, 'BC.pairPermR.csv',row.names = F)
# write.csv(BC.pairPermD, 'BC.pairPermD.csv',row.names = F)
# write.csv(BC.pairPermDR, 'BC.pairPermDR.csv',row.names = F)



## Betadisper of BC (to test the homogeneity of group dispersion) -- Region* ####
# Calculate Bray-Curtis distance matrix
bray_dist <- vegdist(OTMU_bc.cover, method = "bray")
# Test homogeneity of dispersion across regions
dispersion <- betadisper(bray_dist, Factor$Region)

# Perform ANOVA to test for significant differences in dispersion
anova_dispersion <- anova(dispersion)
print(anova_dispersion)
# write.csv(anova_dispersion, 'anova_dispersion_BC_Region.csv',row.names = T)
# Perform permutation test to confirm the ANOVA result
permutation_dispersion <- permutest(dispersion, pairwise = TRUE, permutations = 999)
print(permutation_dispersion)
permutation_dispersion_df <- as.data.frame(permutation_dispersion$tab)
pairwise_results <- as.data.frame(permutation_dispersion$pairwise) # Extract pairwise comparisons
# write.csv(permutation_dispersion_df, 'permutation_dispersion_BC_Region.csv',row.names = T)
# write.csv(pairwise_results, "permutation_dispersion_pairwise_BC_Region.csv", row.names = T)

## Tukey's Honest Significant Differences
(dispersion.HSD <- TukeyHSD(dispersion))
plot(dispersion.HSD)
dispersion.HSD_df <- as.data.frame(dispersion.HSD$group)
# write.csv(dispersion.HSD_df, 'Tukey.HSD_dispersion_BC_Region.csv',row.names = T)

# Visualize dispersion
plot(dispersion)
boxplot(dispersion, main = "Dispersion by Region", ylab = "Distance to Centroid")
# Check distances to centroid for each group
Dispersion <- print(tapply(dispersion$distances, Factor$Region, mean))
Disper_Distance <- as.data.frame(Dispersion)
# write.csv(Disper_Distance, 'Disper_Distance_BC_Region.csv',row.names = T)



## Betadisper of BC (to test the homogeneity of group dispersion) -- Depth ####
# Calculate Bray-Curtis distance matrix
bray_dist <- vegdist(OTMU_bc.cover, method = "bray")
# Test homogeneity of dispersion across depths
dispersion <- betadisper(bray_dist, Factor$Depth)

# Perform ANOVA to test for significant differences in dispersion
anova_dispersion <- anova(dispersion)
print(anova_dispersion)
# write.csv(anova_dispersion, 'anova_dispersion_BC_Depth.csv',row.names = T)
# Perform permutation test to confirm the ANOVA result
permutation_dispersion <- permutest(dispersion, pairwise = TRUE, permutations = 999)
print(permutation_dispersion)

## Tukey's Honest Significant Differences
(dispersion.HSD <- TukeyHSD(dispersion))
plot(dispersion.HSD)
# BC Dispersion in Depth is non-significant!

# Visualize dispersion
plot(dispersion)
boxplot(dispersion, main = "Dispersion by Depth", ylab = "Distance to Centroid")
# Check distances to centroid for each group
Dispersion <- print(tapply(dispersion$distances, Factor$Depth, mean))
Disper_Distance <- as.data.frame(Dispersion)
# write.csv(Disper_Distance, 'Disper_Distance_BC_Depth.csv',row.names = T)



## Betadisper of BC (to test the homogeneity of group dispersion) -- RegionDepth ####
# Calculate Bray-Curtis distance matrix
bray_dist <- vegdist(OTMU_bc.cover, method = "bray")
# Test homogeneity of dispersion across region*depth
dispersion <- betadisper(bray_dist, Factor$RegionDepth)

# Perform ANOVA to test for significant differences in dispersion
anova_dispersion <- anova(dispersion)
print(anova_dispersion)
# write.csv(anova_dispersion, 'anova_dispersion_BC_RegionDepth.csv',row.names = T)
# Perform permutation test to confirm the ANOVA result
permutation_dispersion <- permutest(dispersion, pairwise = TRUE, permutations = 999)
print(permutation_dispersion)

## Tukey's Honest Significant Differences
(dispersion.HSD <- TukeyHSD(dispersion))
plot(dispersion.HSD)
# BC Dispersion in RegionDepth is non-significant!

# Visualize dispersion
plot(dispersion)
boxplot(dispersion, main = "Dispersion by RegionDepth", ylab = "Distance to Centroid")
# Check distances to centroid for each group
Dispersion <- print(tapply(dispersion$distances, Factor$RegionDepth, mean))
Disper_Distance <- as.data.frame(Dispersion)
# write.csv(Disper_Distance, 'Disper_Distance_BC_RegionDepth.csv',row.names = T)

##############################################################################################################################



## Coral Assemblage NMDS---in OTMU level ####
Coral.only <- BC.only %>%
  filter(MajorCategory =='Black Corals'| MajorCategory =='Gorgonian Corals'| MajorCategory =='Soft Corals'| MajorCategory =='Stony Corals')

OTMU_coral.sum <- Coral.only %>% 
  group_by(Region, Location, Depth, MajorCategory, OTMUs) %>% 
  summarise(Number = length(OTMUs))

OTMU_coral.sum$Region <- factor(OTMU_coral.sum$Region, 
                               levels = c('North', 'Green Island','Xiaoliuqiu'))

simi.coral <- OTMU_coral.sum %>%
  mutate(Region = factor(Region, levels = unique(Region)), 
         Location = factor(Location, levels = unique(Location))) %>%
  ungroup() %>%
  group_by(Location,OTMUs) %>% summarize(Abundance = sum(Number))%>%
  pivot_wider(names_from = c(OTMUs), values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  column_to_rownames(var = "Location")

rowSums(simi.coral)
simi.coral <- as.matrix(simi.coral)
OTMU_coral.cover <- proportions(simi.coral, 1)
head(OTMU_coral.cover)
rowSums(OTMU_coral.cover)
OTMU_coral.cover <- as.data.frame(OTMU_coral.cover)

Factor2 <- OTMU_coral.sum
#create factor2 table
Location<-row.names(OTMU_coral.cover)
head(Location)
RegionDepth<-c("North15","North30","North5","North15","North30","North5","North15","North30","North5", 
               "Green Island15","Green Island30","Green Island5","Green Island15","Green Island30","Green Island5","Green Island15","Green Island30","Green Island5", 
               "Xiaoliuqiu15","Xiaoliuqiu30","Xiaoliuqiu5","Xiaoliuqiu15","Xiaoliuqiu30","Xiaoliuqiu5","Xiaoliuqiu15","Xiaoliuqiu30","Xiaoliuqiu5")
Region<-c("North","North","North","North","North","North","North","North","North", 
          "Green Island","Green Island","Green Island","Green Island","Green Island","Green Island","Green Island","Green Island","Green Island", 
          "Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu","Xiaoliuqiu")
Depth<-c("15","30","5","15","30","5","15","30","5",
         "15","30","5","15","30","5","15","30","5",
         "15","30","5","15","30","5","15","30","5")
Factor2<-as.data.frame(cbind(Location,RegionDepth,Region, Depth))
head(Factor2)
Factor2$RegionDepth<-as.factor(Factor2$RegionDepth)
Factor2$Location<-as.factor(Factor2$Location)
Factor2$Region<-as.factor(Factor2$Region)
Factor2$Depth<-as.factor(Factor2$Depth)
#transform cover data due to many 0
OTMU_coral.cover_hell<-decostand(OTMU_coral.cover,"hellinger")

#Meta NMDS
nmds2 <- metaMDS(OTMU_coral.cover, "bray", type='n')

## *indicator species analysis of CA ####
siv2 <- indval(OTMU_coral.cover_hell, Factor2$Region)
gr2 <- siv2$maxcls[siv2$pval <= 0.05]
iv2 <- siv2$indcls[siv2$pval <= 0.05]
pv2 <- siv2$pval[siv2$pval <= 0.05]
fidg2 <- data.frame(group = gr2, indval = iv2, pvalue = pv2)
fidg2 <- fidg2[order(fidg2$group, -fidg2$indval), ]
# write.csv(fidg2, 'indicator species of CA.csv',row.names = T)

# Add group names to fidg2
group_levels2 <- levels(Factor2$Region) # make sure the group fit to right Region (check levels)
fidg2$group_name <- factor(fidg2$group, labels = group_levels2)

# Extract species scores for NMDS plot
species_scores2 <- as.data.frame(vegan::scores(nmds2, display = "species"))
species_scores2$species <- rownames(species_scores2)

# Plot NMDS of CA ####
plot(nmds2, display = "sites", type = 'n', xlab = "NMDS1", ylab = "NMDS2")
# Add hulls for different regions
ordihull(nmds2, groups = Factor2$Region, col = c("#ffc857", "#98c1d9", "#adc178"), lwd = 2)
# Add points for different regions and depths
points(nmds2, "sites", col = "#C2F1F2", pch = 17, select = Factor2$RegionDepth=="North5")
points(nmds2, "sites", col = "#98c1d9", pch = 15, select = Factor2$RegionDepth=="North15")
points(nmds2, "sites", col = "#274c77", pch = 16, select = Factor2$RegionDepth=="North30")
points(nmds2, "sites", col = "#FFEE8C", pch = 17, select = Factor2$RegionDepth=="Green Island5")
points(nmds2, "sites", col = "#ffc857", pch = 15, select = Factor2$RegionDepth=="Green Island15")
points(nmds2, "sites", col = "#ff9f1c", pch = 16, select = Factor2$RegionDepth=="Green Island30")
points(nmds2, "sites", col = "#EAF291", pch = 17, select = Factor2$RegionDepth=="Xiaoliuqiu5")
points(nmds2, "sites", col = "#adc178", pch = 15, select = Factor2$RegionDepth=="Xiaoliuqiu15")
points(nmds2, "sites", col = "#709775", pch = 16, select = Factor2$RegionDepth=="Xiaoliuqiu30")

## *stress of Coral NMDS => nmds2$stress ####
# Add stress value
mtext(paste("Stress:", round(nmds2$stress, 2)), side = 1, line = 3, cex = 0.8, adj = 1)

# wider position
for (i in 1:nrow(fidg2)) {
  species <- rownames(fidg2)[i]
  if (species %in% rownames(species_scores2)) {
    # Adjust x and y coordinates slightly if they are too close to other points (manual approach)
    x <- species_scores2[species, 1]
    y <- species_scores2[species, 2]
    # Example of simple adjustment
    text(x + runif(1, -0.3, 0.3), y + runif(1, -0.5, 0.5), labels = species, cex = 0.5, pos = 4, col = "#979dac")
  }
}
# ***can't add legend :( => i want two legend: Depth symbol + Region color ####


## PERMANOVA of CA ####
# PERMANOVA -> test difference between Region? between Depth? or between Region*Depth?

#Region
PermRegion<-adonis2(OTMU_coral.cover~Region, data=Factor2, method="bray", permutations=999)
pairwise.perm.manova(vegdist(OTMU_coral.cover,"bray"), Factor2$Region, nperm=9999)

#Depth
PermDepth<-adonis2(OTMU_coral.cover~Depth, data=Factor2, method="bray", permutations=999)
pairwise.perm.manova(vegdist(OTMU_coral.cover,"bray"), Factor2$Depth, nperm=9999)

#Depth*Region (interaction between D & R)
Factor2$DepthRegion <- interaction(Factor2$Depth, Factor2$Region)
PermDepthRegion<-adonis2(OTMU_coral.cover~Depth*Region, data=Factor2, method="bray", permutations=999)
pairwise.perm.manova(vegdist(OTMU_coral.cover,"bray"), Factor2$DepthRegion, nperm=9999)


CA.Perm.R <- PermRegion
CA.Perm.D <- PermDepth
CA.Perm.DR <- PermDepthRegion
# write.csv(CA.Perm.R, 'CA.Perm.R.csv',row.names = T)
# write.csv(CA.Perm.D, 'CA.Perm.D.csv',row.names = T)
# write.csv(CA.Perm.DR, 'CA.Perm.DR.csv',row.names = T)

## pairwise PERMANOVA of CA ####
CA.pairPermR <- pairwise.adonis(vegdist(OTMU_coral.cover, "bray"), factors = Factor2$Region, perm = 9999)
CA.pairPermD <- pairwise.adonis(vegdist(OTMU_coral.cover, "bray"), factors = Factor2$Depth, perm = 9999)
CA.pairPermDR <- pairwise.adonis(vegdist(OTMU_coral.cover, "bray"), factors = Factor2$DepthRegion, perm = 9999)
# write.csv(CA.pairPermR, 'CA.pairPermR.csv',row.names = F)
# write.csv(CA.pairPermD, 'CA.pairPermD.csv',row.names = F)
# write.csv(CA.pairPermDR, 'CA.pairPermDR.csv',row.names = F)



## Betadisper of CA (to test the homogeneity of group dispersion) -- Region ####
# Calculate Bray-Curtis distance matrix
bray_dist <- vegdist(OTMU_coral.cover, method = "bray")
# Test homogeneity of dispersion across regions
dispersion <- betadisper(bray_dist, Factor2$Region)

# Perform ANOVA to test for significant differences in dispersion
anova_dispersion <- anova(dispersion)
print(anova_dispersion)
# write.csv(anova_dispersion, 'anova_dispersion_CA_Region.csv',row.names = T)
# Perform permutation test to confirm the ANOVA result
permutation_dispersion <- permutest(dispersion, pairwise = TRUE, permutations = 999)
print(permutation_dispersion)

## Tukey's Honest Significant Differences
(dispersion.HSD <- TukeyHSD(dispersion))
plot(dispersion.HSD)
# CA Dispersion in Region is non-significant!

# Visualize dispersion (optional)
plot(dispersion)
boxplot(dispersion, main = "Dispersion by Region", ylab = "Distance to Centroid")
# Check distances to centroid for each group (optional)
Dispersion <- print(tapply(dispersion$distances, Factor2$Region, mean))
Disper_Distance <- as.data.frame(Dispersion)
# write.csv(Disper_Distance, 'Disper_Distance_CA_Region.csv',row.names = T)



## Betadisper of CA (to test the homogeneity of group dispersion) -- Depth* ####
# Calculate Bray-Curtis distance matrix
bray_dist <- vegdist(OTMU_coral.cover, method = "bray")
# Test homogeneity of dispersion across depths
dispersion <- betadisper(bray_dist, Factor2$Depth)

# Perform ANOVA to test for significant differences in dispersion
anova_dispersion <- anova(dispersion)
print(anova_dispersion)
# write.csv(anova_dispersion, 'anova_dispersion_CA_Depth.csv',row.names = T)
# Perform permutation test to confirm the ANOVA result
permutation_dispersion <- permutest(dispersion, pairwise = TRUE, permutations = 999)
print(permutation_dispersion)
permutation_dispersion_df <- as.data.frame(permutation_dispersion$tab)
pairwise_results <- as.data.frame(permutation_dispersion$pairwise) # Extract pairwise comparisons
# write.csv(permutation_dispersion_df, 'permutation_dispersion_CA_Depth.csv',row.names = T)
# write.csv(pairwise_results, "permutation_dispersion_pairwise_CA_Depth.csv", row.names = T)

## Tukey's Honest Significant Differences
(dispersion.HSD <- TukeyHSD(dispersion))
plot(dispersion.HSD)
dispersion.HSD_df <- as.data.frame(dispersion.HSD$group)
# write.csv(dispersion.HSD_df, 'Tukey.HSD_dispersion_CA_Depth.csv',row.names = T)

# Visualize dispersion (optional)
plot(dispersion)
boxplot(dispersion, main = "Dispersion by Depth", ylab = "Distance to Centroid")
# Check distances to centroid for each group (optional)
Dispersion <- print(tapply(dispersion$distances, Factor2$Depth, mean))
Disper_Distance <- as.data.frame(Dispersion)
# write.csv(Disper_Distance, 'Disper_Distance_CA_Depth.csv',row.names = T)



## Betadisper of CA (to test the homogeneity of group dispersion) -- RegionDepth ####
# Calculate Bray-Curtis distance matrix
bray_dist <- vegdist(OTMU_coral.cover, method = "bray")
# Test homogeneity of dispersion across region*depth
dispersion <- betadisper(bray_dist, Factor2$RegionDepth)

# Perform ANOVA to test for significant differences in dispersion
anova_dispersion <- anova(dispersion)
print(anova_dispersion)
# write.csv(anova_dispersion, 'anova_dispersion_CA_RegionDepth.csv',row.names = T)
# Perform permutation test to confirm the ANOVA result
permutation_dispersion <- permutest(dispersion, pairwise = TRUE, permutations = 999)
print(permutation_dispersion)

## Tukey's Honest Significant Differences
(dispersion.HSD <- TukeyHSD(dispersion))
plot(dispersion.HSD)
# CA Dispersion in RegionDepth is non-significant!

# Visualize dispersion (optional)
plot(dispersion)
boxplot(dispersion, main = "Dispersion by RegionDepth", ylab = "Distance to Centroid")
# Check distances to centroid for each group (optional)
Dispersion <- print(tapply(dispersion$distances, Factor2$RegionDepth, mean))
Disper_Distance <- as.data.frame(Dispersion)
# write.csv(Disper_Distance, 'Disper_Distance_CA_RegionDepth.csv',row.names = T)



### output the new create csv for Beta Diversity analysis ####
# write.csv(OTMU_bc.sum, 'OTMU_bc.sum.csv',row.names = F)
# write.csv(OTMU_coral.sum, 'OTMU_coral.sum.csv',row.names = F)



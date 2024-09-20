### NMDS / PERMANOVA of BC & CA ####
## beta quantitative => NMDS + bray curtis + permanova + simper (indicative species)
library(dplyr); library(tidyr); library(tibble)
library(vegan); library(labdsv); library(RVAideMemoire)

rm (list = ls())
getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R")

### Benthic Community NMDS---in OTU level ####
BC.only <- read.csv("BC.only.csv", header=T, sep=",")

OTU_bc.sum <- BC.only %>% 
  group_by(Region, Location, Depth, MajorCategory, OTUs) %>% 
  summarise(Number = length(OTUs))

OTU_bc.sum$Region <- factor(OTU_bc.sum$Region, 
                            levels = c('North', 'Green Island','Xiaoliuqiu'))

simi.bc <- OTU_bc.sum %>%
  mutate(Region = factor(Region, levels = unique(Region)), 
         Location = factor(Location, levels = unique(Location))) %>%
  ungroup() %>%
  group_by(Location,OTUs) %>% summarize(Abundance = sum(Number))%>%
  pivot_wider(names_from = c(OTUs), values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  column_to_rownames(var = "Location")

rowSums(simi.bc)
simi.bc <- as.matrix(simi.bc)
OTU_bc.cover <- proportions(simi.bc, 1)
head(OTU_bc.cover)
rowSums(OTU_bc.cover)
OTU_bc.cover <- as.data.frame(OTU_bc.cover)

#create Factor table
Factor <- OTU_bc.sum
Location<-row.names(OTU_bc.cover)
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
OTU_bc.cover_hell<-decostand(OTU_bc.cover,"hellinger")

#Meta NMDS
nmds1 <- metaMDS(OTU_bc.cover, "bray", type='n')

## *indicator species analysis of BC ####
siv <- indval(OTU_bc.cover_hell, Factor$Region)
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

# Plot NMDS with adjusted range
plot(nmds1, display = "sites", type = 'n', xlab = "NMDS1", ylab = "NMDS2", xlim = xlim_range, ylim = ylim_range)
# Add points for different regions and depths
points(nmds1, "sites", col = "#e0fbfc", pch = 17, select = Factor$RegionDepth=="North5")
points(nmds1, "sites", col = "#98c1d9", pch = 15, select = Factor$RegionDepth=="North15")
points(nmds1, "sites", col = "#274c77", pch = 16, select = Factor$RegionDepth=="North30")
points(nmds1, "sites", col = "#fff3b0", pch = 17, select = Factor$RegionDepth=="Green Island5")
points(nmds1, "sites", col = "#ffc857", pch = 15, select = Factor$RegionDepth=="Green Island15")
points(nmds1, "sites", col = "#ff9f1c", pch = 16, select = Factor$RegionDepth=="Green Island30")
points(nmds1, "sites", col = "#ecf39e", pch = 17, select = Factor$RegionDepth=="Xiaoliuqiu5")
points(nmds1, "sites", col = "#adc178", pch = 15, select = Factor$RegionDepth=="Xiaoliuqiu15")
points(nmds1, "sites", col = "#709775", pch = 16, select = Factor$RegionDepth=="Xiaoliuqiu30")
# Add hulls for different regions
ordihull(nmds1, groups = Factor$Region, col = c("#ffc857", "#98c1d9", "#adc178"))

## *stress of BC NMDS => nmds1$stress ####
# Add stress value
mtext(paste("Stress:", round(nmds1$stress, 4)), side = 1, line = 3, cex = 0.8, adj = 1)

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

#RegionDepth (need to repeat for Region and Depth!!!)
PermRegionDepth<-adonis(OTU_bc.cover~RegionDepth, data=Factor, method="bray", permutations=999)
PermRegionDepth$aov.tab
pairwise.perm.manova(vegdist(OTU_bc.cover,"bray"), Factor$RegionDepth, nperm=9999)

# BC.Perm.D <- PermDepth$aov.tab
# BC.Perm.R <- PermRegion$aov.tab
# BC.Perm.RD <- PermRegionDepth$aov.tab
# write.csv(BC.Perm.D, 'BC.Perm.D.csv',row.names = F)



## Coral Assemblage NMDS---in OTU level ####
Coral.only <- BC.only %>%
  filter(MajorCategory =='Antipatharian'| MajorCategory =='Gorgonian'| MajorCategory =='Octocoral'| MajorCategory =='Scleractinian')

OTU_coral.sum <- Coral.only %>% 
  group_by(Region, Location, Depth, MajorCategory, OTUs) %>% 
  summarise(Number = length(OTUs))

OTU_coral.sum$Region <- factor(OTU_coral.sum$Region, 
                               levels = c('North', 'Green Island','Xiaoliuqiu'))

simi.coral <- OTU_coral.sum %>%
  mutate(Region = factor(Region, levels = unique(Region)), 
         Location = factor(Location, levels = unique(Location))) %>%
  ungroup() %>%
  group_by(Location,OTUs) %>% summarize(Abundance = sum(Number))%>%
  pivot_wider(names_from = c(OTUs), values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  column_to_rownames(var = "Location")

rowSums(simi.coral)
simi.coral <- as.matrix(simi.coral)
OTU_coral.cover <- proportions(simi.coral, 1)
head(OTU_coral.cover)
rowSums(OTU_coral.cover)
OTU_coral.cover <- as.data.frame(OTU_coral.cover)

Factor2 <- OTU_coral.sum
#create factor2 table
Location<-row.names(OTU_coral.cover)
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
OTU_coral.cover_hell<-decostand(OTU_coral.cover,"hellinger")

#Meta NMDS
nmds2 <- metaMDS(OTU_coral.cover, "bray", type='n')

## *indicator species analysis of CA ####
siv2 <- indval(OTU_coral.cover_hell, Factor2$Region)
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

# Plot NMDS with adjusted range
plot(nmds2, display = "sites", type = 'n', xlab = "NMDS1", ylab = "NMDS2")
# Add points for different regions and depths
points(nmds2, "sites", col = "#e0fbfc", pch = 17, select = Factor2$RegionDepth=="North5")
points(nmds2, "sites", col = "#98c1d9", pch = 15, select = Factor2$RegionDepth=="North15")
points(nmds2, "sites", col = "#274c77", pch = 16, select = Factor2$RegionDepth=="North30")
points(nmds2, "sites", col = "#fff3b0", pch = 17, select = Factor2$RegionDepth=="Green Island5")
points(nmds2, "sites", col = "#ffc857", pch = 15, select = Factor2$RegionDepth=="Green Island15")
points(nmds2, "sites", col = "#ff9f1c", pch = 16, select = Factor2$RegionDepth=="Green Island30")
points(nmds2, "sites", col = "#ecf39e", pch = 17, select = Factor2$RegionDepth=="Xiaoliuqiu5")
points(nmds2, "sites", col = "#adc178", pch = 15, select = Factor2$RegionDepth=="Xiaoliuqiu15")
points(nmds2, "sites", col = "#709775", pch = 16, select = Factor2$RegionDepth=="Xiaoliuqiu30")
# Add hulls for different regions
ordihull(nmds2, groups = Factor2$Region, col = c("#ffc857", "#98c1d9", "#adc178"))

## *stress of Coral NMDS => nmds2$stress ####
# Add stress value
mtext(paste("Stress:", round(nmds2$stress, 4)), side = 1, line = 3, cex = 0.8, adj = 1)

# Add significant indicator species to NMDS plot
for (i in 1:nrow(fidg2)) {
  species <- rownames(fidg2)[i]
  if (species %in% rownames(species_scores2)) {
    text(species_scores2[species, 1], species_scores2[species, 2], labels = species, cex = 0.5, pos = 4, col = "#979dac")
  }
}
# ***can't add legend :( => i want two legend: Depth symbol + Region color ####


## PERMANOVA of CA ####
# PERMANOVA -> test difference between Region? between Depth? or between Region*Depth?

#RegionDepth (need to repeat for Region and Depth!!!)
PermRegionDepth<-adonis(OTU_coral.cover~RegionDepth, data=Factor2, method="bray", permutations=999)
PermRegionDepth$aov.tab
pairwise.perm.manova(vegdist(OTU_coral.cover,"bray"), Factor2$RegionDepth, nperm=9999)

# CA.Perm.D <- PermDepth$aov.tab
# CA.Perm.R <- PermRegion$aov.tab
# CA.Perm.RD <- PermRegionDepth$aov.tab
# write.csv(CA.Perm.D, 'CA.Perm.D.csv',row.names = F)



### output the new create csv for Beta Diversity analysis ####
# write.csv(OTU_bc.sum, 'OTU_bc.sum.csv',row.names = F)
# write.csv(OTU_coral.sum, 'OTU_coral.sum.csv',row.names = F)

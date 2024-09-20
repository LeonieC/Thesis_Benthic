library(ggplot2); library(ggcorrplot); library(readxl); library(corrr); 
library(tidyr); library(dplyr); library(ggbiplot); library(permute); 
library(lattice); library(vegan); library(tidyverse); library(betapart); 
library(rstatix); library(conflicted); library(ggpubr); library(stringr); library(scales)

# install.packages("ggpubr")
conflict_prefer("select", 'dplyr')
conflict_prefer("filter", "dplyr")

rm (list = ls())
getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R")

### Modify Dataset ####
# in another R Script



### Light Attenuation ####
Light.data <- read.csv("PAR data.csv", header=T, sep=",")

Light.data$Region <- factor(Light.data$Region, levels = c('North', 'Green Island','Xiaoliuqiu'))
Light.data$Depth <- factor(Light.data$Depth, levels = c('5', '15','30'))

ggplot(data=Light.data, aes(x=Depth, y=Adjustment, group=Region)) +
  geom_line(aes(color=Region), size=0.8) +
  geom_point(size=1.2) +
  geom_text(aes(label = paste0(sprintf("%.0f", Adjustment * 100), "%")), vjust = -0.8, size = 2.2) +
  scale_color_manual(values=c('#98c1d9','#ffc857','#adc178')) +
  scale_y_continuous(name = "Adjustment value",
                     limits = c(0, 0.35),
                     breaks = seq(0, 0.35, 0.05),
                     labels = percent) +
  labs(title = "Light Attenuation") +
  theme(legend.position="right", 
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 10), legend.text = element_text(size = 8))



### Bayesian Model with Coral Cover---in Genus level ####
# in another R Script



### Benthic Composition---in Major Category level ####
BC.fulldata <- read.csv("BC.fulldata.csv", header=T, sep=",")

MC_sum.numb <- BC.fulldata %>% group_by(Region,Location,Depth,MajorCategory) %>% 
  summarise(Number=length(MajorCategory)) %>% ungroup()

MC_sum.numb$Region <- factor(MC_sum.numb$Region, levels = c('North', 'Green Island','Xiaoliuqiu'))
MC_sum.numb$Depth <- factor(MC_sum.numb$Depth, levels = c('5', '15','30'))

## show total Major Category in the column
# MC_list <- unique(MC_sum.numb$MajorCategory)
# print(MC_list)

MC_sum.numb$MajorCategory <- factor(MC_sum.numb$MajorCategory, 
                                    levels = c('Turf Algae', 'Macroalgae', 'Crustose Coralline Algae',
                                               'Antipatharian', 'Gorgonian', 'Octocoral', 'Scleractinian', 
                                               'Sponge', 'Other Sessile Invertebrates', 'Other Mobile Invertebrates', 'Other Life',
                                               'Unstable Substrate', 'Stable Substrate', 'Shadow', 'Other'))

ggplot(data = MC_sum.numb, mapping = aes(x = Depth, y =Number,fill = MajorCategory)) +
  geom_bar(position = "fill",stat='identity') +
  labs(fill = "Major Category") +
  ylab("Percent") +
  scale_fill_manual(values=c("#dde5b6","#adc178","#c9cba3",
                             "#fee8c8","#fdd49e","#fdbb84","#fc8d59",
                             "#a9d6e5","#61a5c2","#2a6f97","#3d5a80",
                             "#ced4da","#adb5bd","#6c757d","#495057"))+
  facet_wrap(~Region)



### Coral Cover---in Major Category level ####
MC_coral <- MC_sum.numb %>%
  mutate(MajorCategory = ifelse(MajorCategory == 'Antipatharian', 'Antipatharian', 
                                 ifelse(MajorCategory == 'Gorgonian', 'Gorgonian',
                                        ifelse(MajorCategory == 'Octocoral', 'Octocoral',
                                               ifelse(MajorCategory == 'Scleractinian', 'Scleractinian', 'Other'))))) %>%
  group_by(Region, Location, Depth, MajorCategory)%>%
  summarise(Number = sum(Number))

MC_coral.percent <- MC_coral %>%
  group_by(Location) %>% 
  reframe(Total = sum(Number)) %>%
  right_join(., MC_coral) %>%
  group_by(Region, Location, Depth, MajorCategory) %>%
  reframe(Percent = Number/Total)

## output the new create csv
# write.csv(MC_coral.percent, 'MC_coral.percent.csv',row.names = F)

MC_coral.percent$MajorCategory <- factor(MC_coral.percent$MajorCategory, 
                                          levels = c('Other', 'Antipatharian', 'Gorgonian', 'Octocoral', 'Scleractinian'))
MC_coral.percent$Region <- factor(MC_coral.percent$Region, 
                                  levels = c('North', 'Green Island','Xiaoliuqiu'))
MC_coral.percent$Depth <- factor(MC_coral.percent$Depth, 
                                 levels = c('5', '15','30'))

ggplot(data = MC_coral.percent, mapping = aes(x = Depth, y =Percent,fill = MajorCategory)) +
  geom_bar(position = "fill",stat='identity') + 
  labs(fill = "Major Category") +
  scale_fill_manual(values=c('#f8f9fa',"#fff7bc","#fee391","#fec44f","#fe9929")) +
  scale_y_continuous(name = "Coral cover") +
  facet_wrap(~Region)



### Location & Dominant Coral Species---in Genus level ####
BC.only <- read.csv("BC.only.csv", header=T, sep=",")

Coral.only.sum <- BC.only %>% 
  group_by(Region, Location, Depth, MajorCategory, Genus) %>% 
  summarise(Number = length(Genus)) %>%
  filter(MajorCategory =='Antipatharian'| MajorCategory =='Gorgonian'| MajorCategory =='Octocoral'| MajorCategory =='Scleractinian')

Coral.only.sum$Region <- factor(Coral.only.sum$Region, 
                            levels = c('North', 'Green Island','Xiaoliuqiu'))

Coral.only.percent <- Coral.only.sum %>% group_by(Region, Depth) %>% 
  mutate(Percent = Number/sum(Number))

C.domin5.spec <- Coral.only.percent %>% slice_max(Percent, n = 5)  # highest 5 dominant

C.domin5.spec$Region <- factor(C.domin5.spec$Region, 
                               levels = c('North', 'Green Island','Xiaoliuqiu'))
C.domin5.spec$Depth <- factor(C.domin5.spec$Depth, 
                              levels = c('5', '15','30'))
C.domin5.spec$Genus <- factor(C.domin5.spec$Genus, 
                            levels = c('Antipatharia',
                                       'Scleronephthya','Xenia/Heteroxenia','Conglomeratusclera','Sinularia/Cladiella/Klyxum',
                                       'Isopora','Anacropora','Montipora','Porites/Montipora','Porites','Pachyseris',
                                       'Turbinaria','Tubastraea','Tubastraea/Cladopsammia','Pocillopora','Echinophyllia/Oxypora',
                                       'Merulinidae.spp1','Merulinidae.spp2','Cyphastrea','Echinopora','Mycedium','Platygyra'))

ggplot(data = C.domin5.spec, mapping = aes(x = Depth, y =Percent, fill = Genus)) +
  geom_bar(position = "fill", stat='identity') + 
  scale_fill_manual(values=c("#979dac",
                             "#ffd0e9","#ff9ebb","#ff7aa2","#e05780",
                             "#ccd5ae","#abdda4","#66c2a5","#9dd9d2","#92c5de",
                             "#d4a373","#72ddf7","#0096c7","#023e8a","#dab6fc","#9381ff",
                             "#fff7bc","#fee396","#c4ca43","#dec44f","#fe9929","#f47d43")) +
  labs(fill = "Coral Genus") +
  theme(legend.position="right",
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 12), axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 10), axis.text.y = element_text(size = 10),
        legend.title = element_text(size = 10), legend.text = element_text(size = 8)) +
  facet_wrap(~Region)



### Alpha Diversity (OTU richness)---in OTU level => per Transect OTUs ####
BC.only <- read.csv("BC.only.csv", header=T, sep=",")

OTU_sum.numb <- BC.only %>% 
  mutate(Number = 1) %>% 
  group_by(Region, Site, Transect, Depth, MajorCategory, Genus, OTUs) %>% 
  summarise(Number=sum(Number))

OTU_sum.otu <- OTU_sum.numb %>% mutate(OTU = 1) %>% 
  group_by(Region,Site,Transect,Depth) %>%
  summarise(OTU=sum(OTU)) %>%
  ungroup()

OTU_sum.otu$Region <- factor(OTU_sum.otu$Region, 
                             levels = c('North', 'Green Island','Xiaoliuqiu'))
OTU_sum.otu$Depth <- factor(OTU_sum.otu$Depth, 
                            levels = c('5', '15','30'))

ggplot(OTU_sum.otu, aes(x = Depth, y = OTU, fill = Depth))+
  geom_boxplot(outliers = F)+
  geom_jitter(size=0.8)+
  facet_wrap(~Region)

## *significant test of OTU ####
# all region
shapiro.test(OTU_sum.otu$OTU) #normality test-> p<0.05 means not Normal Distribution
kru_OTU <- kruskal_test(OTU ~ Region, data = OTU_sum.otu) #non parametric test 無母數(n<30)
wil_OTU <- wilcox_test(OTU ~ Region, p.adjust.method = 'bonferroni', data = OTU_sum.otu) %>% #pairwise wilcoxon test
  add_xy_position(x = "Region") 

ggplot(OTU_sum.otu, aes(x=Region, y=OTU))+
  geom_boxplot()+
  stat_pvalue_manual(wil_OTU, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)

#North
OTU_sum.otu_N <- OTU_sum.otu %>% filter(Region == 'North')
shapiro.test(OTU_sum.otu_N$OTU) 
kru_OTU_N <- kruskal_test(OTU ~ Depth, data = OTU_sum.otu_N) 
wil_OTU_N <- wilcox_test(OTU ~ Depth, p.adjust.method = 'bonferroni', data = OTU_sum.otu_N) %>% 
  add_xy_position(x = "Depth")

ggplot(OTU_sum.otu_N, aes(x=Depth, y=OTU))+
  geom_boxplot()+
  stat_pvalue_manual(wil_OTU_N, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)+
  ggtitle("North")+
  theme(plot.title = element_text(hjust = 0.5))

#GI
OTU_sum.otu_GI <- OTU_sum.otu %>% filter(Region == 'Green Island')
shapiro.test(OTU_sum.otu_GI$OTU) 
kru_OTU_GI <- kruskal_test(OTU ~ Depth, data = OTU_sum.otu_GI) 
wil_OTU_GI <- wilcox_test(OTU ~ Depth, p.adjust.method = 'bonferroni', data = OTU_sum.otu_GI) %>% 
  add_xy_position(x = "Depth")

ggplot(OTU_sum.otu_GI, aes(x=Depth, y=OTU))+
  geom_boxplot()+
  stat_pvalue_manual(wil_OTU_GI, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)+
  ggtitle("Green Island")+
  theme(plot.title = element_text(hjust = 0.5))

#XLQ
OTU_sum.otu_XLQ <- OTU_sum.otu %>% filter(Region == 'Xiaoliuqiu')
shapiro.test(OTU_sum.otu_XLQ$OTU) 
kru_OTU_XLQ <- kruskal_test(OTU ~ Depth, data = OTU_sum.otu_XLQ) 
wil_OTU_XLQ <- wilcox_test(OTU ~ Depth, p.adjust.method = 'bonferroni', data = OTU_sum.otu_XLQ) %>% 
  add_xy_position(x = "Depth")

ggplot(OTU_sum.otu_XLQ, aes(x=Depth, y=OTU))+
  geom_boxplot()+
  stat_pvalue_manual(wil_OTU_XLQ, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)+
  ggtitle("Xiaoliuqiu")+
  theme(plot.title = element_text(hjust = 0.5))

## *OTU SD ####
OTU_sd.otu <- OTU_sum.otu %>% 
  group_by(Region,Depth) %>% 
  summarise(SD = sd(OTU, na.rm = F),
            Mean=mean(OTU))

# combine the OTU richness plot with significant annotation in each region & depth
# add Region column
wil_OTU_N$Region <- 'North'
wil_OTU_GI$Region <- 'Green Island'
wil_OTU_XLQ$Region <- 'Xiaoliuqiu'

# combine into a single data frame
combined_wil <- bind_rows(wil_OTU_N, wil_OTU_GI, wil_OTU_XLQ)
combined_wil$Region <- factor(combined_wil$Region, 
                             levels = c('North', 'Green Island','Xiaoliuqiu'))

ggplot(OTU_sum.otu, aes(x = Depth, y = OTU)) +
  geom_boxplot(outliers = FALSE) +
  geom_jitter(size = 0.8) +
  facet_wrap(~Region) +
  stat_pvalue_manual(combined_wil, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
  theme(plot.title = element_text(hjust = 0.5))

# ***change Depth color of significant annotation plot ####
ggboxplot(OTU_sum.otu, x = "Depth", y = "OTU", fill = "Depth",
          palette = c("#F8AF97", "#F5BD7F", "#FFEA91"), facet.by = "Region",
          outlier.shape = NA) +
  geom_jitter(size = 0.8) +
  stat_pvalue_manual(combined_wil, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, step.increase = 0.01) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"))  # Adjust the size of the legend keys (symbols)

# ***add background for significant annotation plot ####
# just for try... "#F8AF97", "#F5BD7F", "#FFEA91" / "#CAD2C5", "#84A98C", "#52796F"
ggboxplot(OTU_sum.otu, x = "Depth", y = "OTU", fill = "Depth",
          palette = c("#E7B800", "#608c38", "#557b9d"), facet.by = "Region",
          outlier.shape = NA) +
  geom_jitter(size = 0.8) +
  stat_pvalue_manual(combined_wil, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  theme(legend.position = "right",
        legend.text = element_text(size = 8),
        legend.title = element_text(size = 10),
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_line(color = "lightgray", linewidth = 0.5), # Major grid lines
        panel.grid.minor = element_line(color = "lightgray", linewidth = 0.5)) # Minor grid lines



## beta quantitative => NMDS + bray curtis + permanova + simper (indicative species)
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

Factor <- OTU_bc.sum #create factor table
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
nmds1<-metaMDS(OTU_bc.cover,"bray",type='n')
plot(nmds1, display="sites", type='n')
points(nmds1, "sites",col="#e0fbfc",pch=17,select=Factor$RegionDepth=="North5")
points(nmds1, "sites",col="#98c1d9",pch=15,select=Factor$RegionDepth=="North15")
points(nmds1, "sites",col="#274c77",pch=16,select=Factor$RegionDepth=="North30")
points(nmds1, "sites",col="#fff3b0",pch=17,select=Factor$RegionDepth=="Green Island5")
points(nmds1, "sites",col="#ffc857",pch=15,select=Factor$RegionDepth=="Green Island15")
points(nmds1, "sites",col="#ff9f1c",pch=16,select=Factor$RegionDepth=="Green Island30")
points(nmds1, "sites",col="#ecf39e",pch=17,select=Factor$RegionDepth=="Xiaoliuqiu5")
points(nmds1, "sites",col="#adc178",pch=15,select=Factor$RegionDepth=="Xiaoliuqiu15")
points(nmds1, "sites",col="#709775",pch=16,select=Factor$RegionDepth=="Xiaoliuqiu30")
ordihull(nmds1,groups=Factor$Region, col=c("#ffc857","#98c1d9","#adc178"))
## *stress of BC NMDS => nmds1$stress ####
mtext(paste("Stress:", round(nmds1$stress, 4)), side = 1, line = 3, cex = 0.8, adj = 1)
# text(nmds1, "species", cex=0.5)
# ***can't add legend :( => i want two legend: Depth symbol + Region color ####

## PERMANOVA of BC ####
# PERMANOVA -> test difference between Region? between Depth? or between Region*Depth?
library(RVAideMemoire)

#RegionDepth (need to repeat for Region and Depth!!!)
PermRegionDepth<-adonis(OTU_bc.cover~RegionDepth, data=Factor, method="bray", permutations=999)
PermRegionDepth$aov.tab
pairwise.perm.manova(vegdist(OTU_bc.cover,"bray"), Factor$RegionDepth, nperm=9999)

## *legend indicative species of BC ####
library(labdsv)
siv<-indval(OTU_bc.cover_hell, Factor$Region)
gr<-siv$maxcls[siv$pval<=0.05]
iv<-siv$indcls[siv$pval<=0.05]
pv<-siv$pval[siv$pval<=0.05]
fidg<-data.frame(group=gr, indval=iv, pvalue=pv)
(fidg<-fidg[order(fidg$group,-fidg$indval),]) # ()is for print the result
# ***how to add indicative species with group? what is the group?? ####



### try3 ####
library(vegan)
library(labdsv)

# Perform NMDS
nmds1 <- metaMDS(OTU_bc.cover, "bray", type='n')

# Perform indicator species analysis
siv <- indval(OTU_bc.cover_hell, Factor$Region)
gr <- siv$maxcls[siv$pval <= 0.05]
iv <- siv$indcls[siv$pval <= 0.05]
pv <- siv$pval[siv$pval <= 0.05]
fidg <- data.frame(group = gr, indval = iv, pvalue = pv)
fidg <- fidg[order(fidg$group, -fidg$indval), ]

# Add group names to fidg
group_levels <- levels(Factor$Region)
fidg$group_name <- factor(fidg$group, labels = group_levels)

print(fidg)

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
points(nmds1, "sites", col = "#e0fbfc", pch = 17, select = Factor$Region == "North" & Factor$Depth == 5)
points(nmds1, "sites", col = "#98c1d9", pch = 15, select = Factor$Region == "North" & Factor$Depth == 15)
points(nmds1, "sites", col = "#274c77", pch = 16, select = Factor$Region == "North" & Factor$Depth == 30)
points(nmds1, "sites", col = "#fff3b0", pch = 17, select = Factor$Region == "Green Island" & Factor$Depth == 5)
points(nmds1, "sites", col = "#ffc857", pch = 15, select = Factor$Region == "Green Island" & Factor$Depth == 15)
points(nmds1, "sites", col = "#ff9f1c", pch = 16, select = Factor$Region == "Green Island" & Factor$Depth == 30)
points(nmds1, "sites", col = "#ecf39e", pch = 17, select = Factor$Region == "Xiaoliuqiu" & Factor$Depth == 5)
points(nmds1, "sites", col = "#adc178", pch = 15, select = Factor$Region == "Xiaoliuqiu" & Factor$Depth == 15)
points(nmds1, "sites", col = "#709775", pch = 16, select = Factor$Region == "Xiaoliuqiu" & Factor$Depth == 30)

# Add hulls for different regions
ordihull(nmds1, groups = Factor$Region, col = c("#ffc857", "#98c1d9", "#adc178"))

# Add stress value
mtext(paste("Stress:", round(nmds1$stress, 4)), side = 1, line = 3, cex = 0.8, adj = 1)

# Add significant indicator species to NMDS plot
for (i in 1:nrow(fidg)) {
  species <- rownames(fidg)[i]
  if (species %in% rownames(species_scores)) {
    text(species_scores[species, 1], species_scores[species, 2], labels = species, cex = 0.5, pos = 4, col = "#979dac")
  }
}
### try3 ####



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

#Meta nMDS
nmds2<-metaMDS(OTU_coral.cover,"bray",type='n')
plot(nmds2, display="sites", type='n')
points(nmds2, "sites",col="#e0fbfc",pch=17,select=Factor2$RegionDepth=="North5")
points(nmds2, "sites",col="#98c1d9",pch=15,select=Factor2$RegionDepth=="North15")
points(nmds2, "sites",col="#274c77",pch=16,select=Factor2$RegionDepth=="North30")
points(nmds2, "sites",col="#fff3b0",pch=17,select=Factor2$RegionDepth=="Green Island5")
points(nmds2, "sites",col="#ffc857",pch=15,select=Factor2$RegionDepth=="Green Island15")
points(nmds2, "sites",col="#ff9f1c",pch=16,select=Factor2$RegionDepth=="Green Island30")
points(nmds2, "sites",col="#ecf39e",pch=17,select=Factor2$RegionDepth=="Xiaoliuqiu5")
points(nmds2, "sites",col="#adc178",pch=15,select=Factor2$RegionDepth=="Xiaoliuqiu15")
points(nmds2, "sites",col="#709775",pch=16,select=Factor2$RegionDepth=="Xiaoliuqiu30")
ordihull(nmds2,groups=Factor2$Region, col=c("#ffc857","#98c1d9","#adc178"))
## *stress of Coral NMDS => nmds2$stress ####
mtext(paste("Stress:", round(nmds2$stress, 4)), side = 1, line = 3, cex = 0.8, adj = 1)
# text(nmds2, "species", cex=0.5)
# ***can't add legend :( => i want two legend: Depth symbol + Region color ####

## PERMANOVA of Coral ####
# PERMANOVA -> test difference between Region? between Depth? or between Region*Depth?
library(RVAideMemoire)

#RegionDepth (need to repeat for Region and Depth!!!)
PermRegionDepth<-adonis(OTU_coral.cover~RegionDepth, data=Factor2, method="bray", permutations=999)
PermRegionDepth$aov.tab
pairwise.perm.manova(vegdist(OTU_coral.cover,"bray"), Factor2$RegionDepth, nperm=9999)

## *legend indicative species of Coral ####
library(labdsv)
siv2<-indval(OTU_coral.cover_hell, Factor2$Region)
gr2<-siv2$maxcls[siv2$pval<=0.05]
iv2<-siv2$indcls[siv2$pval<=0.05]
pv2<-siv2$pval[siv2$pval<=0.05]
fidg2<-data.frame(group=gr2, indval=iv2, pvalue=pv2)
(fidg2<-fidg2[order(fidg2$group,-fidg2$indval),])
# ***how to add indicative species with group? what is the group?? ####



### Beta Diversity: Turnover & Nestedness => heatmap---in OTU level ####
### *all Benthic Community data ####
beta.BC <- OTU_bc.sum %>%
  mutate(Reg.Dep = paste(Region,Depth,sep = "_")) %>%
  mutate(Region = factor(Region, levels = unique(Region)), 
         Reg.Dep = factor(Reg.Dep, levels = unique(Reg.Dep))) %>%
  ungroup() %>%
  group_by(Reg.Dep,OTUs) %>% summarize(Abundance = sum(Number)) %>%
  pivot_wider(names_from = c(OTUs), values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  column_to_rownames(var = "Reg.Dep")

OTU<-colnames(beta.BC)
location <- row.names(beta.BC)
beta.BC <- matrix(unlist(beta.BC),nrow = 9, ncol=230)
rownames(beta.BC) <- location 
colnames(beta.BC) <-OTU
beta_jar_mat <- ifelse(beta.BC,beta.BC[] >1,0)

beta.multi(beta_jar_mat, index.family ="jaccard")
beta_location <- beta.pair(beta_jar_mat, index.family="jaccard")

## Turnover of BC ####
beta_location_tu <- beta_location$beta.jtu
dst <- data.matrix(beta_location_tu)

dst_tu_nor <- dst[1:9,1:9]  # 1:3->N / 4:6->GI / 7:9->XLQ
dim_tu_nor <- ncol(dst_tu_nor)

x_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")

image(1:dim_tu_nor, 1:dim_tu_nor, dst_tu_nor, axes = F, xlab="", ylab="")
axis(1, at = 1:dim_tu_nor, labels = x_labels, cex.axis = 1.0, font.axis = 2, las=1)
axis(2, at = 1:dim_tu_nor, labels = y_labels, cex.axis = 1.0, font.axis = 2, las=1)
text(expand.grid(1:dim_tu_nor, 1:dim_tu_nor), sprintf("%0.2f", dst_tu_nor), cex=1.5, font=2)

## Nestedness of BC ####
beta_location_ne <- beta_location$beta.jne
dst <- data.matrix(beta_location_ne)

dst_ne_nor <- dst[1:9,1:9] # 1:3->N / 4:6->GI / 7:9->XLQ
dim_ne_nor <- ncol(dst_ne_nor)

x_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")

image(1:dim_ne_nor, 1:dim_ne_nor, dst_ne_nor, axes = F, xlab="", ylab="")
axis(1, at = 1:dim_ne_nor, labels = x_labels, cex.axis = 1.0, font.axis = 2, las=1)
axis(2, at = 1:dim_ne_nor, labels = y_labels, cex.axis = 1.0, font.axis = 2, las=1)
text(expand.grid(1:dim_ne_nor, 1:dim_ne_nor), sprintf("%0.2f", dst_ne_nor), cex=1.5, font=2)

## Beta Diversity of BC ####
beta_location_ac <- beta_location$beta.jac
dst <- data.matrix(beta_location_ac)

dst_ac_nor <- dst[1:9,1:9] # 1:3->N / 4:6->GI / 7:9->XLQ
dim_ac_nor <- ncol(dst_ac_nor)

x_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")

image(1:dim_ac_nor, 1:dim_ac_nor, dst_ac_nor, axes = F, xlab="", ylab="")
axis(1, at = 1:dim_ac_nor, labels = x_labels, cex.axis = 1.0, font.axis = 2, las=1)
axis(2, at = 1:dim_ac_nor, labels = y_labels, cex.axis = 1.0, font.axis = 2, las=1)
text(expand.grid(1:dim_ac_nor, 1:dim_ac_nor), sprintf("%0.2f", dst_ac_nor), cex=1.5, font=2)



### *only Coral Assemblage data ####
beta.CA <- OTU_coral.sum %>%
  mutate(Reg.Dep = paste(Region,Depth,sep = "_")) %>%
  mutate(Region = factor(Region, levels = unique(Region)), 
         Reg.Dep = factor(Reg.Dep, levels = unique(Reg.Dep))) %>%
  ungroup() %>%
  group_by(Reg.Dep,OTUs) %>% summarize(Abundance = sum(Number)) %>%
  pivot_wider(names_from = c(OTUs), values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  column_to_rownames(var = "Reg.Dep")

OTU<-colnames(beta.CA)
location <- row.names(beta.CA)
beta.CA <- matrix(unlist(beta.CA),nrow = 9, ncol=144)
rownames(beta.CA) <- location 
colnames(beta.CA) <-OTU
beta_jar_mat.ca <- ifelse(beta.CA,beta.CA[] >1,0)

beta.multi(beta_jar_mat.ca, index.family ="jaccard")
beta_location.ca <- beta.pair(beta_jar_mat.ca, index.family="jaccard")

## Turnover of CA ####
beta_location.ca_tu <- beta_location.ca$beta.jtu
dst <- data.matrix(beta_location.ca_tu)

dst_tu_nor <- dst[1:9,1:9]  # 1:3->N / 4:6->GI / 7:9->XLQ
dim_tu_nor <- ncol(dst_tu_nor)

x_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")

image(1:dim_tu_nor, 1:dim_tu_nor, dst_tu_nor, axes = F, xlab="", ylab="")
axis(1, at = 1:dim_tu_nor, labels = x_labels, cex.axis = 1.0, font.axis = 2, las=1)
axis(2, at = 1:dim_tu_nor, labels = y_labels, cex.axis = 1.0, font.axis = 2, las=1)
text(expand.grid(1:dim_tu_nor, 1:dim_tu_nor), sprintf("%0.2f", dst_tu_nor), cex=1.5, font=2)

## Nestedness of CA ####
beta_location.ca_ne <- beta_location.ca$beta.jne
dst <- data.matrix(beta_location.ca_ne)

dst_ne_nor <- dst[1:9,1:9] # 1:3->N / 4:6->GI / 7:9->XLQ
dim_ne_nor <- ncol(dst_ne_nor)

x_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")

image(1:dim_ne_nor, 1:dim_ne_nor, dst_ne_nor, axes = F, xlab="", ylab="")
axis(1, at = 1:dim_ne_nor, labels = x_labels, cex.axis = 1.0, font.axis = 2, las=1)
axis(2, at = 1:dim_ne_nor, labels = y_labels, cex.axis = 1.0, font.axis = 2, las=1)
text(expand.grid(1:dim_ne_nor, 1:dim_ne_nor), sprintf("%0.2f", dst_ne_nor), cex=1.5, font=2)

## Beta Diversity of CA ####
beta_location.ca_ac <- beta_location.ca$beta.jac
dst <- data.matrix(beta_location.ca_ac)

dst_ac_nor <- dst[1:9,1:9] # 1:3->N / 4:6->GI / 7:9->XLQ
dim_ac_nor <- ncol(dst_ac_nor)

x_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")

image(1:dim_ac_nor, 1:dim_ac_nor, dst_ac_nor, axes = F, xlab="", ylab="")
axis(1, at = 1:dim_ac_nor, labels = x_labels, cex.axis = 1.0, font.axis = 2, las=1)
axis(2, at = 1:dim_ac_nor, labels = y_labels, cex.axis = 1.0, font.axis = 2, las=1)
text(expand.grid(1:dim_ac_nor, 1:dim_ac_nor), sprintf("%0.2f", dst_ac_nor), cex=1.5, font=2)

#######

#######


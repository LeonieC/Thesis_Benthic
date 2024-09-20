library(ggplot2); library(ggcorrplot); library(readxl); library(corrr); 
library(tidyr); library(dplyr); library(ggbiplot); library(permute); 
library(lattice); library(vegan); library(tidyverse); library(betapart); 
library(rstatix); library(conflicted); library(ggpubr); library(stringr); library(scales)

# install.packages("ggpubr")
conflict_prefer("select", 'dplyr')
conflict_prefer("filter", "dplyr")

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



### Bayesian Model with Coral Cover ####
# in another R Script



### Benthic Composition in Major Category level ####
BC.fulldata <- read.csv("BC.fulldata.csv", header=T, sep=",")

MC_sum.numb <- BC.fulldata %>% group_by(Region,Location,Depth,MajorCategory) %>% 
  summarise(Number=length(MajorCategory)) %>% ungroup()

MC_sum.numb$Region <- factor(MC_sum.numb$Region, levels = c('North', 'Green Island','Xiaoliuqiu'))
MC_sum.numb$Depth <- factor(MC_sum.numb$Depth, levels = c('5', '15','30'))

# show total Major Category in the column
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

#######

#######

# rename FG_sum.numb to MC_sum.numb

### coral and algae change position
FG_sum.numb$Taxonomy.Group <- factor(FG_sum.numb$Taxonomy.Group, 
                                     levels = c('Antipatharia', 'Gorgonian', 'Octocoral', 'Scleractinian',
                                                'Turf Algae', 'Macroalgae', 'Crustose Coralline Algae',
                                                'Sponge', 'Other Sessile Benthos', 'Other Life',
                                                'Soft Substrate', 'Hard Substrate', 'Other'))

ggplot(data = FG_sum.numb, mapping = aes(x = Region, y =Number,fill = Taxonomy.Group)) +
  geom_bar(position = "fill",stat='identity') + 
  scale_fill_manual(values=c("#fee8c8","#fdd49e","#fdbb84","#fc8d59",
                             "#dde5b6","#adc178","#c9cba3",
                             "#a9d6e5","#61a5c2","#2a6f97",
                             "#ced4da","#adb5bd","#6c757d"))+
  facet_wrap(~Depth)

ggplot(data = FG_sum.numb, mapping = aes(x = Region, y =Number,fill = Taxonomy.Group)) +
  geom_bar(position = "fill",stat='identity') + 
  scale_fill_manual(values=c("#fff7bc","#fee391","#fec44f","#fe9929",
                             "#dde5b6","#adc178","#c9cba3",
                             "#a9d6e5","#61a5c2","#2a6f97",
                             "#ced4da","#adb5bd","#6c757d"))+
  facet_wrap(~Depth)

##################################################################################

## Coral Cover ####
FG_coral <- FG_sum.numb %>%
  mutate(Taxonomy.Group = ifelse(Taxonomy.Group == 'Antipatharia', 'Antipatharia', 
                                 ifelse(Taxonomy.Group == 'Gorgonian', 'Gorgonian',
                                        ifelse(Taxonomy.Group == 'Octocoral', 'Octocoral',
                                               ifelse(Taxonomy.Group == 'Scleractinian', 'Scleractinian', 'Other'))))) %>%
  group_by(Location, Taxonomy.Group, Region, Depth)%>%
  summarise(Number = sum(Number))

FG_coral.percent <- FG_coral %>%
  group_by(Location) %>% 
  reframe(Total = sum(Number)) %>%
  right_join(., FG_coral) %>%
  group_by(Location, Taxonomy.Group, Region, Depth) %>%
  reframe(Percent = Number/Total) 


## output the new create csv
# write.csv(FG_coral.percent, 'FG_coral.percent.csv',row.names = F)

FG_coral.percent$Taxonomy.Group <- factor(FG_coral.percent$Taxonomy.Group, 
                                          levels = c('Other', 'Antipatharia', 'Gorgonian', 'Octocoral', 'Scleractinian'))
FG_coral.percent$Region <- factor(FG_coral.percent$Region, 
                                  levels = c('North', 'Green Island','Xiaoliuqiu'))
FG_coral.percent$Depth <- factor(FG_coral.percent$Depth, 
                                 levels = c('5', '15','30'))

ggplot(data = FG_coral.percent, mapping = aes(x = Depth, y =Percent,fill = Taxonomy.Group)) +
  geom_bar(position = "fill",stat='identity') + 
  scale_fill_manual(values=c('#f8f9fa',"#fff7bc","#fee391","#fec44f","#fe9929"))+
  facet_wrap(~Region)


#####???
FG_Region <- FG_coral.percent %>% ungroup() %>%
  group_by(Taxonomy.Group, Region, Depth)

ggplot(data = FG_Region, mapping = aes(x = Region, y =Percent,fill = Taxonomy.Group)) +
  geom_bar(position = "fill",stat='identity') + 
  scale_fill_manual(values=c('#f8f9fa',"#fff7bc","#fee391","#fec44f","#fe9929"))+
  facet_wrap(~Depth) + 
  geom_text(aes(label = paste(round(Percent,digits = 2)*100,'%')))

##################################################################################

# tidyr, dplyr --> combine two file for species name
BC_code <- read.csv("short_code.csv", header=T, sep=",") %>%
  group_by(Label, Genus) %>%
  summarise(C= sum(Cover))

BC_genus.name <- BC.data %>% 
  left_join(., BC_code, by = 'Label') %>%
  select(-C)

# the data only contain the living things!!
BC_life.genus.name <- BC_genus.name %>%
  filter(Taxonomy.Group !='Hard Substrate' & Taxonomy.Group !='Soft Substrate')


# Location & Dominant Species ####
c_spec.only <- BC_life.genus.name %>% 
  group_by(Region, Location, Depth, Taxonomy.Group, Genus) %>% 
  summarise(Number = length(Genus)) %>%
  filter(Taxonomy.Group =='Antipatharia'| Taxonomy.Group =='Gorgonian'| Taxonomy.Group =='Octocoral'| Taxonomy.Group =='Scleractinian')

c_spec.only$Region <- factor(c_spec.only$Region, 
                             levels = c('North', 'Green Island','Xiaoliuqiu'))

c_spec.only.percent <- c_spec.only %>% group_by(Region, Depth) %>% 
  mutate(ratio = Number/sum(Number))

c_domin5.spec <- c_spec.only.percent %>% slice_max(ratio, n = 5)  # highest 5 dominant


c_domin5.spec$Region <- factor(c_domin5.spec$Region, 
                               levels = c('North', 'Green Island','Xiaoliuqiu'))
c_domin5.spec$Depth <- factor(c_domin5.spec$Depth, 
                              levels = c('5', '15','30'))

# choose those species bigger than 10% #
# mutate(Dominant = ifelse(ratio > 0.1, Genus, "Other"))
# check sum is 100% or not #
# summarise(sum(ratio))


ggplot(data = c_domin5.spec, mapping = aes(x = Depth, y =ratio,fill = Genus)) +
  geom_bar(position = "fill", stat='identity') + 
  scale_fill_manual(values=c("#ccd5ae","#d53e4f","#f46d43","#92c5de","#38b000",
                             "#c1fba4","#e6f598","#abdda4","#66c2a5","#fe9929",
                             "#d4a373","#9970ab","#fb9a99","#de77ae","#b79ced",
                             '#979dac',"#fff7bc","#fee391","#fec44f","#4393c3",
                             "#9d4edd","#709775","#2166ac"))+
  facet_wrap(~Region)

# color ####
# "#f46d43","#d53e4f","#fb9a99","#de77ae","#b79ced","#9970ab","#9d4edd"
# "#a4c379","#40916c","#709795","#c1fba4","#e6f598","#709775","#38b000","#A5CD9F"
# "#fee391","#f4cd43","#fec44f"
# "#4393c3","#2166ac","#caf0f8","#90e0ef"

## benthic dominant ####
bcl <- BC_life.genus.name %>% 
  group_by(Region, Location, Depth, Label, Taxonomy.Group, Genus) %>% 
  summarise(Number = length(Label)) 

bcl$Region <- factor(bcl$Region, levels = c('North', 'Green Island','Xiaoliuqiu'))

bcl.percent <- bcl %>% group_by(Region, Depth) %>% 
  mutate(ratio = Number/sum(Number))

bcl_domin10.spec <- bcl.percent %>% slice_max(ratio, n = 10)  # highest 10 dominant


bcl_domin10.spec$Region <- factor(bcl_domin10.spec$Region, 
                                  levels = c('North', 'Green Island','Xiaoliuqiu'))
bcl_domin10.spec$Depth <- factor(bcl_domin10.spec$Depth, 
                                 levels = c('5', '15','30'))

ggplot(data = bcl_domin10.spec, mapping = aes(x = Depth, y =ratio,fill = Label)) +
  geom_bar(position = "fill", stat='identity') + 
  scale_fill_manual(values=c("#f46d43","#d53e4f","#e0afa0","#92c5de","#fee340",
                             "#fc8d59","#e6f598","#abdda4","#fee360","#fe9929",
                             "#fb9a99","#9970ab","#fee391","#f15bb5","#b79ced",
                             '#979dac',"#fff7bc","#d4a373","#fec44f","#c1fba4",
                             "#ccd5ae","#709775","#2166ac"))+
  facet_wrap(~Region)


##################################################################################
### OLD ONE ###
## Region & Depth similarity => NMDS ##

### Coral Assemblage
simi_C.only <- c_spec.only %>%
  mutate(Region = factor(Region, levels = unique(Region)), 
         Location = factor(Location, levels = unique(Location))) %>%
  ungroup() %>%
  group_by(Location,Genus) %>% summarize(Abundance = sum(Number))%>%
  pivot_wider(names_from = c(Genus), values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  column_to_rownames(var = "Location")

Abund.nmds <- metaMDS(simi_C.only, k = 2, distance = 'bray')

png(filename = "NMDS_Coral Assemblage.png", width = 6, height = 4, units = "in", res = 300)
par(mar = c(4, 4, 2, 2))
ordiplot(Abund.nmds, type = "none")
ordihull(Abund.nmds, groups = as.factor(rep(c('North', 'GI', 'XLQ'), each = 9)), draw = 'polygon', 
         col = c("#9d4edd","#2166ac","#709775"))

points(Abund.nmds, display = "sites", col = rep(c("#2166ac","#9d4edd","#709775"), each = 9),
       pch = rep(c(15,16,17), times = 9), cex = 0.75)

legend('topright', legend = c('North', 'GI', 'XLQ'),
       text.col = c("#2166ac","#9d4edd","#709775"),
       bty = 'n', title = 'Region',title.adj = 0.75 ,title.col = "black", cex = 0.75)
legend('topleft', title = 'Depth',title.col = "black",legend = c('5 m', '15 m', '30 m'), pch = c(17,15,16),
       bty = 'n', cex = 0.75)
dev.off()


### Benthic Community
c.spec_sum <- BC_life.genus.name %>% 
  group_by(Region, Location, Depth, Taxonomy.Group, Genus) %>% 
  summarise(Number = length(Genus))

c.spec_sum$Region <- factor(c.spec_sum$Region, 
                            levels = c('North', 'Green Island','Xiaoliuqiu'))

simi_BC <- c.spec_sum %>%
  mutate(Region = factor(Region, levels = unique(Region)), 
         Location = factor(Location, levels = unique(Location))) %>%
  ungroup() %>%
  group_by(Location,Genus) %>% summarize(Abundance = sum(Number))%>%
  pivot_wider(names_from = c(Genus), values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  column_to_rownames(var = "Location")

BC.nmds <- metaMDS(simi_BC, k = 2, distance = 'bray')

png(filename = "NMDS_Benthic Community.png", width = 6, height = 4, units = "in", res = 300)
par(mar = c(4, 4, 2, 2))
ordiplot(BC.nmds, type = "none")
ordihull(BC.nmds, groups = as.factor(rep(c('North', 'GI', 'XLQ'), each = 9)), draw = 'polygon', 
         col = c("#9d4edd","#2166ac","#709775"))

points(BC.nmds, display = "sites", col = rep(c("#2166ac","#9d4edd","#709775"), each = 9),
       pch = rep(c(15,16,17), times = 9), cex = 0.75)

legend('topright', legend = c('North', 'GI', 'XLQ'),
       text.col = c("#2166ac","#9d4edd","#709775"),
       bty = 'n', title = 'Region',title.adj = 0.75 ,title.col = "black", cex = 0.75)
legend('topleft', title = 'Depth',title.col = "black",legend = c('5 m', '15 m', '30 m'), pch = c(17,15,16),
       bty = 'n', cex = 0.75)
dev.off()


##################################################################################
### NEW ONE ###
## alpha diversity of OTU ?? (OTU richness) ####
BC_life.genus.name
#OLD--Log strightly on Number
OTU_sum.numb <- BC_life.genus.name %>% 
  mutate(Number = 1) %>% 
  group_by(Region,Transect,Depth,Taxonomy.Group,Genus,Label) %>% 
  summarise(Number=sum(Number)) %>%
  mutate(Log = log10(Number)) %>%
  filter(Taxonomy.Group !='Hard Substrate' & Taxonomy.Group !='Soft Substrate') %>%
  ungroup()
OTU_sum.numb$Region <- factor(OTU_sum.numb$Region, 
                              levels = c('North', 'Green Island','Xiaoliuqiu'))
OTU_sum.numb$Depth <- factor(OTU_sum.numb$Depth, 
                             levels = c('5', '15','30'))
ggplot(OTU_sum.numb, aes(x = Region, y = Log, fill = Depth))+
  geom_boxplot()
## OLD ##

##NEW--Log on y axis
OTU_sum.numb <- BC_life.genus.name %>% 
  mutate(Number = 1) %>% 
  group_by(Region,Transect,Depth,Taxonomy.Group,Genus,Label) %>% 
  summarise(Number=sum(Number)) %>%
  filter(Taxonomy.Group !='Hard Substrate' & Taxonomy.Group !='Soft Substrate') %>%
  ungroup()
OTU_sum.numb$Region <- factor(OTU_sum.numb$Region, 
                              levels = c('North', 'Green Island','Xiaoliuqiu'))
OTU_sum.numb$Depth <- factor(OTU_sum.numb$Depth, 
                             levels = c('5', '15','30'))
ggplot(OTU_sum.numb, aes(x = Region, y = Number, fill = Depth))+
  geom_boxplot()

p <- ggplot(OTU_sum.numb, aes(x = Region, y = Number, fill = Depth))+
  geom_boxplot()
p + scale_y_continuous(trans = "log10")
## NEW ##

# NEW 2 -- Correct one per transect OTUs #
OTU_sum.numb <- BC_life.genus.name %>% 
  mutate(Number = 1) %>% 
  group_by(Region,Site,Transect,Depth,Taxonomy.Group,Genus,Label) %>% 
  summarise(Number=sum(Number)) %>%
  filter(Taxonomy.Group !='Hard Substrate' & Taxonomy.Group !='Soft Substrate')
OTU_sum.otu <- OTU_sum.numb %>% mutate(OTUs = 1) %>% 
  group_by(Region,Site,Transect,Depth) %>%
  summarise(OTUs=sum(OTUs)) %>%
  ungroup()
OTU_sum.otu$Region <- factor(OTU_sum.otu$Region, 
                             levels = c('North', 'Green Island','Xiaoliuqiu'))
OTU_sum.otu$Depth <- factor(OTU_sum.otu$Depth, 
                            levels = c('5', '15','30'))

ggplot(OTU_sum.otu, aes(x = Depth, y = OTUs, fill = Depth))+
  geom_boxplot(outliers = F)+
  geom_jitter(size=0.8)+
  facet_wrap(~Region)
# NEW 2 #

## significant test ####
# all region
shapiro.test(OTU_sum.otu$OTUs) #normality test-> p<0.05 means not Normal Distribution
kru_OTU <- kruskal_test(OTUs ~ Region, data = OTU_sum.otu) #non parametric test 無母數(n<30)
wil_OTU <- wilcox_test(OTUs ~ Region, p.adjust.method = 'bonferroni', data = OTU_sum.otu) %>% #pairwise wilcoxon test
  add_xy_position(x = "Region") 

ggplot(OTU_sum.otu, aes(x=Region, y=OTUs))+
  geom_boxplot()+
  stat_pvalue_manual(wil_OTU, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)


#North
OTU_sum.otu_N <- OTU_sum.otu %>% filter(Region == 'North')
shapiro.test(OTU_sum.otu_N$OTUs) 
kru_OTU_N <- kruskal_test(OTUs ~ Depth, data = OTU_sum.otu_N) 
wil_OTU_N <- wilcox_test(OTUs ~ Depth, p.adjust.method = 'bonferroni', data = OTU_sum.otu_N) %>% 
  add_xy_position(x = "Depth")

ggplot(OTU_sum.otu_N, aes(x=Depth, y=OTUs))+
  geom_boxplot()+
  stat_pvalue_manual(wil_OTU_N, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)+
  ggtitle("North")+
  theme(plot.title = element_text(hjust = 0.5))

#GI
OTU_sum.otu_GI <- OTU_sum.otu %>% filter(Region == 'Green Island')
shapiro.test(OTU_sum.otu_GI$OTUs) 
kru_OTU_GI <- kruskal_test(OTUs ~ Depth, data = OTU_sum.otu_GI) 
wil_OTU_GI <- wilcox_test(OTUs ~ Depth, p.adjust.method = 'bonferroni', data = OTU_sum.otu_GI) %>% 
  add_xy_position(x = "Depth")

ggplot(OTU_sum.otu_GI, aes(x=Depth, y=OTUs))+
  geom_boxplot()+
  stat_pvalue_manual(wil_OTU_GI, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)+
  ggtitle("Green Island")+
  theme(plot.title = element_text(hjust = 0.5))


#XLQ
OTU_sum.otu_XLQ <- OTU_sum.otu %>% filter(Region == 'Xiaoliuqiu')
shapiro.test(OTU_sum.otu_XLQ$OTUs) 
kru_OTU_XLQ <- kruskal_test(OTUs ~ Depth, data = OTU_sum.otu_XLQ) 
wil_OTU_XLQ <- wilcox_test(OTUs ~ Depth, p.adjust.method = 'bonferroni', data = OTU_sum.otu_XLQ) %>% 
  add_xy_position(x = "Depth")

ggplot(OTU_sum.otu_XLQ, aes(x=Depth, y=OTUs))+
  geom_boxplot()+
  stat_pvalue_manual(wil_OTU_XLQ, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)+
  ggtitle("Xiaoliuqiu")+
  theme(plot.title = element_text(hjust = 0.5))


### Overall OTUs ####
OTU_overall.numb <- BC_life.genus.name %>% 
  ungroup()%>%
  mutate(Number = 1) %>% 
  group_by(Label) %>% 
  summarise(Number=sum(Number))

OTU_overall.sum <- OTU_overall.numb %>% mutate(OTUs = 1) %>% 
  group_by(Label) %>%
  summarise(OTUs=sum(OTUs))
OTUs.total<- OTU_overall.sum%>%
  summarise(OTUs=sum(OTUs))

## OTU SD ####
#OLD
OTU_sd.numb <- OTU_sum.numb %>% 
  group_by(Region,Depth) %>% 
  summarise(SD = sd(Number, na.rm = F),
            Mean=mean(Number))

#NEW
OTU_sd.otu <- OTU_sum.otu %>% 
  group_by(Region,Depth) %>% 
  summarise(SD = sd(OTUs, na.rm = F),
            Mean=mean(OTUs))

##################################################################################

### beta quantitative => NMDS + bray curtis + permanova + simper (indicative species)

## Benthic Community (OTU) NMDS ####
rowSums(simi_BC)
simi_BC <- as.matrix(simi_BC)
OTU_cover <- proportions(simi_BC, 1)
head(OTU_cover)
rowSums(OTU_cover)
OTU_cover <- as.data.frame(OTU_cover)
Coral_cover<- #dataset with benthi compositin in % for only CORAL
  
  Factor <- c.spec_sum
#create factor table
Location<-row.names(OTU_cover)
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
OTU_cover_hell<-decostand(OTU_cover,"hellinger")

"#e0fbfc","#98c1d9","#274c77"
"#ecf39e","#adc178","#709775"
"#fff3b0","#ffc857","#ff9f1c"

#Meta nMDS
nmds1<-metaMDS(OTU_cover,"bray",type='n')
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
# nmds1$stress
ordihull(nmds1,groups=Factor$Region, col=c("#ffc857","#98c1d9","#adc178"))
text(nmds1, "species", cex=0.5)

#Permanova-> difference btw region? betwen depth? or between depth*region
#region
library(RVAideMemoire) # Region, Depth, RegionDepth
PermRegionDepth<-adonis(OTU_cover~Depth, data=Factor, method="bray", permutations=999)
PermRegionDepth$aov.tab
pairwise.perm.manova(vegdist(OTU_cover,"bray"), Factor$Depth, nperm=9999)

#legendre indicative species
library(labdsv)
siv<-indval(OTU_cover_hell, Factor$RegionDepth)
gr<-siv$maxcls[siv$pval<=0.05]
iv<-siv$indcls[siv$pval<=0.05]
pv<-siv$pval[siv$pval<=0.05]
fidg<-data.frame(group=gr, indval=iv, pvalue=pv)
(fidg<-fidg[order(fidg$group,-fidg$indval),])


## Coral Assemblage NMDS ####
rowSums(simi_C.only)
simi_C.only <- as.matrix(simi_C.only)
Coral_cover <- proportions(simi_C.only, 1)
head(Coral_cover)
rowSums(Coral_cover)
Coral_cover <- as.data.frame(Coral_cover)

Factor <- c.spec_sum
#create factor table
Location<-row.names(Coral_cover)
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
Coral_cover_hell<-decostand(Coral_cover,"hellinger")

"#e0fbfc","#98c1d9","#274c77"
"#ecf39e","#adc178","#709775"
"#fff3b0","#ffc857","#ff9f1c"

#Meta nMDS
nmds2<-metaMDS(Coral_cover,"bray",type='n')
plot(nmds2, display="sites", type='n')
points(nmds2, "sites",col="#e0fbfc",pch=17,select=Factor$RegionDepth=="North5")
points(nmds2, "sites",col="#98c1d9",pch=15,select=Factor$RegionDepth=="North15")
points(nmds2, "sites",col="#274c77",pch=16,select=Factor$RegionDepth=="North30")
points(nmds2, "sites",col="#fff3b0",pch=17,select=Factor$RegionDepth=="Green Island5")
points(nmds2, "sites",col="#ffc857",pch=15,select=Factor$RegionDepth=="Green Island15")
points(nmds2, "sites",col="#ff9f1c",pch=16,select=Factor$RegionDepth=="Green Island30")
points(nmds2, "sites",col="#ecf39e",pch=17,select=Factor$RegionDepth=="Xiaoliuqiu5")
points(nmds2, "sites",col="#adc178",pch=15,select=Factor$RegionDepth=="Xiaoliuqiu15")
points(nmds2, "sites",col="#709775",pch=16,select=Factor$RegionDepth=="Xiaoliuqiu30")
# nmds2$stress
ordihull(nmds2,groups=Factor$Region, col=c("#ffc857","#98c1d9","#adc178"))
text(nmds2, "species", cex=0.5)

#Permanova-> difference btw region? betwen depth? or between depth*region
#region
library(RVAideMemoire) # Region, Depth, RegionDepth
PermRegionDepth<-adonis(Coral_cover~Depth, data=Factor, method="bray", permutations=999)
PermRegionDepth$aov.tab
pairwise.perm.manova(vegdist(Coral_cover,"bray"), Factor$Depth, nperm=9999)

#legendre indicative species
library(labdsv)
siv<-indval(Coral_cover_hell, Factor$RegionDepth)
gr<-siv$maxcls[siv$pval<=0.05]
iv<-siv$indcls[siv$pval<=0.05]
pv<-siv$pval[siv$pval<=0.05]
fidg<-data.frame(group=gr, indval=iv, pvalue=pv)
(fidg<-fidg[order(fidg$group,-fidg$indval),])

###??
region_group_mapping <- data.frame(Region = unique(ShallowSiteFactor$Region2), Group = gr)
coral.ind <-  multipatt(Coral_cover_hell, Factor$RegionDepth, restcomb = 1, control = how(nperm=999)) # calculate indval values of fish species


##################################################################################

### betapart turnover & nestedness => heatmap ####

# levels=c("North_5","North_15","North_30", "Green Island_5","Green Island_15","Green Island_30", "Xiaoliuqiu_5","Xiaoliuqiu_15","Xiaoliuqiu_30")

### all BC ####
beta_bc <- c.spec_sum %>%
  mutate(Reg.Dep = paste(Region,Depth,sep = "_")) 

beta_BC <- beta_bc %>%
  mutate(Region = factor(Region, levels = unique(Region)), 
         Reg.Dep = factor(Reg.Dep, levels = unique(Reg.Dep))) %>%
  ungroup() %>%
  group_by(Reg.Dep,Genus) %>% summarize(Abundance = sum(Number)) %>%
  pivot_wider(names_from = c(Genus), values_from = Abundance, values_fill = list(Abundance = 0))

write.csv(beta_BC, 'beta_BC.csv',row.names = F)
beta_BC.re <- read.csv("beta_BC.csv", header=T, sep=",") %>%
  column_to_rownames(var = "Reg.Dep")


OTU<-colnames(beta_BC.re)
location <- row.names(beta_BC.re)
beta_BC.re <- matrix(unlist(beta_BC.re),nrow = 9, ncol=130)
rownames(beta_BC.re) <- location 
colnames(beta_BC.re) <-OTU
beta_jar_mat <- ifelse(beta_BC.re,beta_BC.re[] >1,0)



beta.multi(beta_jar_mat, index.family ="jaccard")
beta_location <- beta.pair(beta_jar_mat, index.family="jaccard")

## Turnover ####
beta_location_tu <- beta_location$beta.jtu

dst <- data.matrix(beta_location_tu)
dim <- ncol(dst)

image(1:dim, 1:dim, dst, axes = FALSE, xlab="", ylab="")
axis(1, 1:dim, cex.axis = 0.5, las=1)
axis(2, 1:dim, cex.axis = 0.5, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%0.1f", dst), cex=1)

dst_tu_nor <- dst[1:9,1:9]  # 1:3->N / 4:6->GI / 7:9->XLQ
nor_lab <- colnames(dst_tu_nor) 
dim_tu_nor <- ncol(dst_tu_nor)
image(1:dim_tu_nor, 1:dim_tu_nor, dst_tu_nor, axes = F, xlab="", ylab="")
axis(1, 1:dim_tu_nor, nor_lab, cex.axis = 0.5, las=1)
axis(2, 1:dim_tu_nor, nor_lab, cex.axis = 0.5, las=1)
text(expand.grid(1:dim_tu_nor, 1:dim_tu_nor), sprintf("%0.2f", dst_tu_nor), cex=1)

## Nestedness ####
beta_location_ne <- beta_location$beta.jne

dst <- data.matrix(beta_location_ne)
dim <- ncol(dst)

image(1:dim, 1:dim, dst, axes = FALSE, xlab="", ylab="")
axis(1, 1:dim, cex.axis = 0.5, las=1)
axis(2, 1:dim, cex.axis = 0.5, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%0.1f", dst), cex=0.6)

dst_ne_nor <- dst[1:9,1:9] # 1:3->N / 4:6->GI / 7:9->XLQ
nor_lab <- colnames(dst_ne_nor) 
dim_ne_nor <- ncol(dst_ne_nor)
image(1:dim_ne_nor, 1:dim_ne_nor, dst_ne_nor, axes = F, xlab="", ylab="")
axis(1, 1:dim_ne_nor, nor_lab, cex.axis = 0.5, las=1)
axis(2, 1:dim_ne_nor, nor_lab, cex.axis = 0.5, las=1)
text(expand.grid(1:dim_ne_nor, 1:dim_ne_nor), sprintf("%0.2f", dst_ne_nor), cex=1)

## beta diversity ####
beta_location_ac <- beta_location$beta.jac

dst <- data.matrix(beta_location_ac)
dim <- ncol(dst)

image(1:dim, 1:dim, dst, axes = FALSE, xlab="", ylab="")
axis(1, 1:dim, cex.axis = 0.5, las=1)
axis(2, 1:dim, cex.axis = 0.5, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%0.1f", dst), cex=0.6)

dst_ne_nor <- dst[1:9,1:9] # 1:3->N / 4:6->GI / 7:9->XLQ
nor_lab <- colnames(dst_ne_nor) 
dim_ne_nor <- ncol(dst_ne_nor)
image(1:dim_ne_nor, 1:dim_ne_nor, dst_ne_nor, axes = F, xlab="", ylab="")
axis(1, 1:dim_ne_nor, nor_lab, cex.axis = 0.5, las=1)
axis(2, 1:dim_ne_nor, nor_lab, cex.axis = 0.5, las=1)
text(expand.grid(1:dim_ne_nor, 1:dim_ne_nor), sprintf("%0.2f", dst_ne_nor), cex=1)


### only CA ####
beta_ca <- c_spec.only %>%
  mutate(Reg.Dep = paste(Region,Depth,sep = "_")) 

beta_CA <- beta_ca %>%
  mutate(Region = factor(Region, levels = unique(Region)), 
         Reg.Dep = factor(Reg.Dep, levels = unique(Reg.Dep))) %>%
  ungroup() %>%
  group_by(Reg.Dep,Genus) %>% summarize(Abundance = sum(Number)) %>%
  pivot_wider(names_from = c(Genus), values_from = Abundance, values_fill = list(Abundance = 0))

write.csv(beta_CA, 'beta_CA.csv',row.names = F)
beta_CA.re <- read.csv("beta_CA.csv", header=T, sep=",") %>%
  column_to_rownames(var = "Reg.Dep")


OTU<-colnames(beta_CA.re)
location <- row.names(beta_CA.re)
beta_CA.re <- matrix(unlist(beta_CA.re),nrow = 9, ncol=130)
rownames(beta_CA.re) <- location 
colnames(beta_CA.re) <-OTU
beta_jar_mat.ca <- ifelse(beta_CA.re,beta_CA.re[] >1,0)



beta.multi(beta_jar_mat.ca, index.family ="jaccard")
beta_location.ca <- beta.pair(beta_jar_mat.ca, index.family="jaccard")

## Turnover ####
beta_location.ca_tu <- beta_location.ca$beta.jtu

dst <- data.matrix(beta_location.ca_tu)
dim <- ncol(dst)

image(1:dim, 1:dim, dst, axes = FALSE, xlab="", ylab="")
axis(1, 1:dim, cex.axis = 0.5, las=1)
axis(2, 1:dim, cex.axis = 0.5, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%0.1f", dst), cex=1)

dst_tu_nor <- dst[1:9,1:9]  # 1:3->N / 4:6->GI / 7:9->XLQ
nor_lab <- colnames(dst_tu_nor) 
dim_tu_nor <- ncol(dst_tu_nor)
image(1:dim_tu_nor, 1:dim_tu_nor, dst_tu_nor, axes = F, xlab="", ylab="")
axis(1, 1:dim_tu_nor, nor_lab, cex.axis = 0.5, las=1)
axis(2, 1:dim_tu_nor, nor_lab, cex.axis = 0.5, las=1)
text(expand.grid(1:dim_tu_nor, 1:dim_tu_nor), sprintf("%0.2f", dst_tu_nor), cex=1)

## Nestedness ####
beta_location.ca_ne <- beta_location.ca$beta.jne

dst <- data.matrix(beta_location.ca_ne)
dim <- ncol(dst)

image(1:dim, 1:dim, dst, axes = FALSE, xlab="", ylab="")
axis(1, 1:dim, cex.axis = 0.5, las=1)
axis(2, 1:dim, cex.axis = 0.5, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%0.1f", dst), cex=0.6)

dst_ne_nor <- dst[1:9,1:9] # 1:3->N / 4:6->GI / 7:9->XLQ
nor_lab <- colnames(dst_ne_nor) 
dim_ne_nor <- ncol(dst_ne_nor)
image(1:dim_ne_nor, 1:dim_ne_nor, dst_ne_nor, axes = F, xlab="", ylab="")
axis(1, 1:dim_ne_nor, nor_lab, cex.axis = 0.5, las=1)
axis(2, 1:dim_ne_nor, nor_lab, cex.axis = 0.5, las=1)
text(expand.grid(1:dim_ne_nor, 1:dim_ne_nor), sprintf("%0.2f", dst_ne_nor), cex=1)

## beta diversity ####
beta_location.ca_ac <- beta_location.ca$beta.jac

dst <- data.matrix(beta_location.ca_ac)
dim <- ncol(dst)

image(1:dim, 1:dim, dst, axes = FALSE, xlab="", ylab="")
axis(1, 1:dim, cex.axis = 0.5, las=1)
axis(2, 1:dim, cex.axis = 0.5, las=1)
text(expand.grid(1:dim, 1:dim), sprintf("%0.1f", dst), cex=0.6)

dst_ne_nor <- dst[1:9,1:9] # 1:3->N / 4:6->GI / 7:9->XLQ
nor_lab <- colnames(dst_ne_nor) 
dim_ne_nor <- ncol(dst_ne_nor)
image(1:dim_ne_nor, 1:dim_ne_nor, dst_ne_nor, axes = F, xlab="", ylab="")
axis(1, 1:dim_ne_nor, nor_lab, cex.axis = 0.5, las=1)
axis(2, 1:dim_ne_nor, nor_lab, cex.axis = 0.5, las=1)
text(expand.grid(1:dim_ne_nor, 1:dim_ne_nor), sprintf("%0.2f", dst_ne_nor), cex=1)







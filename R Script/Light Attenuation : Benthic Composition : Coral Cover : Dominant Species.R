### Light Attenuation / Benthic Composition / Coral Cover / Dominant Species ####
library(ggplot2); library(dplyr); library(scales)

rm (list = ls())
getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R")

### Light Attenuation ####
PAR.diving.data <- read.csv("PAR raw data diving.csv", header=T, sep=",")
PAR.surface.data <- read.csv("PAR raw data surface.csv", header=T, sep=",")

## *PAR SD ####
PAR_sd <- PAR.diving.data %>% 
  group_by(Region, Depth) %>% 
  summarise(SD = sd(PAR, na.rm = F),
            Mean=mean(PAR))

## *Attenuation Plot ####
PAR.diving.mean <- PAR.diving.data %>%
  group_by(Region, Site, Depth) %>% 
  summarise(D.Mean=mean(PAR)) %>%
  ungroup()

PAR.surface.mean <- PAR.surface.data %>%
  group_by(Region, Site, Depth) %>% 
  summarise(S.Mean=mean(PAR)) %>%
  ungroup()

PAR.full <- PAR.diving.mean %>%
  left_join(select(PAR.surface.mean, Site, Depth, S.Mean), by = c("Site", "Depth")) %>%
  mutate(Adjustment = D.Mean / S.Mean) %>%
  group_by(Region, Depth) %>%
  summarise(D.Mean = mean(D.Mean),
            S.Mean = mean(S.Mean),
            Adjustment = mean(Adjustment))

PAR.full$Region <- factor(PAR.full$Region, levels = c('North', 'Green Island','Xiaoliuqiu'))
PAR.full$Depth <- factor(PAR.full$Depth, levels = c('5', '15','30'))

ggplot(data=PAR.full, aes(x=Depth, y=Adjustment, group=Region)) +
  geom_line(aes(color=Region), size=0.8) +
  geom_point(size=1.2) +
  geom_text(aes(label = paste0(sprintf("%.1f", Adjustment * 100), "%")), vjust = -0.8, size = 2.2) +
  scale_color_manual(values=c('#98c1d9','#ffc857','#adc178')) +
  scale_y_continuous(name = "Adjustment value",
                     limits = c(0, 0.35),
                     breaks = seq(0, 0.35, 0.05),
                     labels = percent) +
  theme(legend.position="right", 
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title.x = element_text(size = 14), axis.title.y = element_text(size = 14),
        axis.text.x = element_text(size = 12), axis.text.y = element_text(size = 12),
        legend.title = element_text(size = 10), legend.text = element_text(size = 8))

PAR.RS.sd <- PAR.full %>%
  group_by(Region) %>%
  summarise(SD = sd(S.Mean, na.rm = F),
            R.S.Mean = mean(S.Mean))



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


# ***calculate % for BC ####
MC_RD.sum.numb <- MC_sum.numb %>%
  mutate(RegDep = paste0(Region, '_', Depth)) %>%
  group_by(Region, RegDep, Depth ,MajorCategory) %>%
  summarise(Number = sum(Number))
MC_BC.percent <- MC_RD.sum.numb %>%
  group_by(RegDep) %>%
  reframe(Total = sum(Number)) %>%
  right_join(., MC_RD.sum.numb) %>%
  group_by(Region, RegDep, Depth ,MajorCategory) %>%
  summarise(Percent = Number/Total)

MC_BC.percent <- MC_RD.sum.numb %>%
  group_by(RegDep) %>%
  reframe(Total = sum(Number)) %>%
  right_join(MC_RD.sum.numb, by = c("RegDep")) %>%
  mutate(Percent = Number / Total) %>%
  group_by(Region, Depth ,MajorCategory) %>%
  summarise(Percent = Percent,
            Percent2 = sprintf("%.2f%%", Percent * 100), 
            .groups = 'drop')
# write.csv(MC_BC.percent, 'MC_BC.percent.csv',row.names = F)



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

# ***calculate % for coral ####
MC_RD.C.sum.numb <- MC_coral %>%
  mutate(RegDep = paste0(Region, '_', Depth)) %>%
  group_by(Region, RegDep, Depth ,MajorCategory) %>%
  summarise(Number = sum(Number))
MC_C.percent <- MC_RD.C.sum.numb %>%
  group_by(RegDep) %>%
  reframe(Total = sum(Number)) %>%
  right_join(MC_RD.C.sum.numb, by = c("RegDep")) %>%
  mutate(Percent = Number / Total) %>%
  group_by(Region, Depth ,MajorCategory) %>%
  summarise(Percent = Percent,
            Percent2 = sprintf("%.2f%%", Percent * 100), 
            .groups = 'drop')
# write.csv(MC_C.percent, 'MC_C.percent.csv',row.names = F)



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

# write.csv(C.domin5.spec, 'C.5dominant.genus.csv',row.names = F)


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



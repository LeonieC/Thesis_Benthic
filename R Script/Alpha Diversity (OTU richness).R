### Alpha Diversity (OTU richness) ####
library(ggplot2); library(dplyr); library(rstatix); library(ggpubr)

rm (list = ls())
getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R")

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



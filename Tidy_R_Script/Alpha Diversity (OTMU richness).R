### Alpha Diversity (OTMU richness) ####
library(ggplot2); library(dplyr); library(rstatix); library(ggpubr)

rm (list = ls())
getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file")

### Alpha Diversity (OTMU richness)---in OTMU level => per Transect OTMUs ####
BC.only <- read.csv("BC.only.csv", header=T, sep=",")

OTMU_sum.numb <- BC.only %>% 
  mutate(Number = 1) %>% 
  group_by(Region, Site, Transect, Depth, MajorCategory, OTUs, OTMUs) %>% 
  summarise(Number=sum(Number))

OTMU_sum.otmu <- OTMU_sum.numb %>% mutate(OTMU = 1) %>% 
  group_by(Region,Site,Transect,Depth) %>%
  summarise(OTMU=sum(OTMU)) %>%
  ungroup()

OTMU_sum.otmu$Region <- factor(OTMU_sum.otmu$Region, 
                             levels = c('North', 'Green Island','Xiaoliuqiu'))
OTMU_sum.otmu$Depth <- factor(OTMU_sum.otmu$Depth, 
                            levels = c('5', '15','30'))

ggplot(OTMU_sum.otmu, aes(x = Depth, y = OTMU, fill = Depth))+
  geom_boxplot(outliers = F)+
  geom_jitter(size=0.8)+
  facet_wrap(~Region)

## *significant test of OTMU ####
# all region
shapiro.test(OTMU_sum.otmu$OTMU) #normality test-> p<0.05 means not Normal Distribution
kru_OTMU <- kruskal_test(OTMU ~ Region, data = OTMU_sum.otmu) #non parametric test 無母數(n<30)
# kruskal.test(OTMU ~ Region, data = OTMU_sum.otmu) # same as above one
wil_OTMU <- wilcox_test(OTMU ~ Region, p.adjust.method = 'bonferroni', data = OTMU_sum.otmu) %>% #pairwise wilcoxon test
  add_xy_position(x = "Region") %>%
  mutate(p.adj.signif = case_when(p.adj < 0.001 ~ "***", TRUE ~ p.adj.signif))

ggplot(OTMU_sum.otmu, aes(x=Region, y=OTMU))+
  geom_boxplot()+
  stat_pvalue_manual(wil_OTMU, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)

#North
OTMU_sum.otmu_N <- OTMU_sum.otmu %>% filter(Region == 'North')
shapiro.test(OTMU_sum.otmu_N$OTMU) 
kru_OTMU_N <- kruskal_test(OTMU ~ Depth, data = OTMU_sum.otmu_N) 
wil_OTMU_N <- wilcox_test(OTMU ~ Depth, p.adjust.method = 'bonferroni', data = OTMU_sum.otmu_N) %>% 
  add_xy_position(x = "Depth") %>%
  mutate(p.adj.signif = case_when(p.adj < 0.001 ~ "***", TRUE ~ p.adj.signif))

ggplot(OTMU_sum.otmu_N, aes(x=Depth, y=OTMU))+
  geom_boxplot()+
  stat_pvalue_manual(wil_OTMU_N, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)+
  ggtitle("North")+
  theme(plot.title = element_text(hjust = 0.5))

#GI
OTMU_sum.otmu_GI <- OTMU_sum.otmu %>% filter(Region == 'Green Island')
shapiro.test(OTMU_sum.otmu_GI$OTMU) 
kru_OTMU_GI <- kruskal_test(OTMU ~ Depth, data = OTMU_sum.otmu_GI) 
wil_OTMU_GI <- wilcox_test(OTMU ~ Depth, p.adjust.method = 'bonferroni', data = OTMU_sum.otmu_GI) %>% 
  add_xy_position(x = "Depth") %>%
  mutate(p.adj.signif = case_when(p.adj < 0.001 ~ "***", TRUE ~ p.adj.signif))

ggplot(OTMU_sum.otmu_GI, aes(x=Depth, y=OTMU))+
  geom_boxplot()+
  stat_pvalue_manual(wil_OTMU_GI, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)+
  ggtitle("Green Island")+
  theme(plot.title = element_text(hjust = 0.5))

#XLQ
OTMU_sum.otmu_XLQ <- OTMU_sum.otmu %>% filter(Region == 'Xiaoliuqiu')
shapiro.test(OTMU_sum.otmu_XLQ$OTMU) 
kru_OTMU_XLQ <- kruskal_test(OTMU ~ Depth, data = OTMU_sum.otmu_XLQ) 
wil_OTMU_XLQ <- wilcox_test(OTMU ~ Depth, p.adjust.method = 'bonferroni', data = OTMU_sum.otmu_XLQ) %>% 
  add_xy_position(x = "Depth") %>%
  mutate(p.adj.signif = case_when(p.adj < 0.001 ~ "***", TRUE ~ p.adj.signif))

ggplot(OTMU_sum.otmu_XLQ, aes(x=Depth, y=OTMU))+
  geom_boxplot()+
  stat_pvalue_manual(wil_OTMU_XLQ, label = "p.adj.signif", tip.length = 0.01, hide.ns = T)+
  ggtitle("Xiaoliuqiu")+
  theme(plot.title = element_text(hjust = 0.5))


## *OTMU SD ####
OTMU_sd.otmu <- OTMU_sum.otmu %>% 
  group_by(Region,Depth) %>% 
  summarise(SD = sd(OTMU, na.rm = F),
            Mean=mean(OTMU))
# write.csv(OTMU_sd.otmu, 'OTMU richness SD.csv',row.names = T)

# combine the OTMU richness plot with significant annotation in each region & depth
# add Region column
wil_OTMU_N$Region <- 'North'
wil_OTMU_GI$Region <- 'Green Island'
wil_OTMU_XLQ$Region <- 'Xiaoliuqiu'

# combine into a single data frame
combined_wil <- bind_rows(wil_OTMU_N, wil_OTMU_GI, wil_OTMU_XLQ)
combined_wil$Region <- factor(combined_wil$Region, 
                              levels = c('North', 'Green Island','Xiaoliuqiu'))

# ***change Depth color of significant annotation plot ####
OTMU_sum.otmu$Region_Depth <- paste(OTMU_sum.otmu$Region, OTMU_sum.otmu$Depth, sep = "_")
OTMU_sum.otmu <- OTMU_sum.otmu %>%
  mutate(Region_Depth = gsub("^North_", "North Taiwan_", Region_Depth))
OTMU_sum.otmu$Region_Depth <- factor(OTMU_sum.otmu$Region_Depth, 
                                     levels = c('North Taiwan_5', 'North Taiwan_15', 'North Taiwan_30',
                                                'Green Island_5', 'Green Island_15', 'Green Island_30',
                                                'Xiaoliuqiu_5', 'Xiaoliuqiu_15', 'Xiaoliuqiu_30'))

ggboxplot(OTMU_sum.otmu, x = "Depth", y = "OTMU", fill = "Region_Depth",  # Fill by updated Region_Depth
          palette = c("#e0fbfc", "#98c1d9", "#205c81", "#fff3b0", "#ffc857", "#ff9f1c", "#ecf39e", "#adc178", "#709775"), 
          facet.by = "Region",
          outlier.shape = NA) +
  geom_jitter(size = 0.8) +
  stat_pvalue_manual(combined_wil, label = "p.adj.signif", tip.length = 0.01, hide.ns = TRUE, step.increase = 0.01) +
  scale_y_continuous(name = "OTMU",
                     limits = c(10, 50),
                     breaks = seq(10, 50, 10)) +
  xlab("Depth (m)") +
  theme(legend.position = "right",
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        strip.text = element_text(size = 13, face = "bold"),
        legend.text = element_text(size = 9),
        legend.title = element_text(size = 11),
        legend.key.size = unit(0.5, "cm"),
        panel.grid.major = element_line(color = "#e9ecef"),
        panel.grid.minor = element_line(color = "#e9ecef")) +
  facet_wrap(~Region, 
             labeller = labeller(Region = c("North" = "North Taiwan", 
                                            "Green Island" = "Green Island", 
                                            "Xiaoliuqiu" = "Xiaoliuqiu")))  # Modify facet labels


# sum info ####
OTMU_sum <- OTMU_sum.otmu %>%
  group_by(Region_Depth, Depth) %>%
  summarise(OTMU=sum(OTMU)) %>%
  mutate(Transect = case_when(
    Depth == 5 ~ 15,  # 105 pictures in 5 m in each site
    Depth == 15 ~ 15, # 105 pictures in 15 m in each site
    Depth == 30 ~ 9)) %>%
  mutate(aveg_OTMU = OTMU/Transect)



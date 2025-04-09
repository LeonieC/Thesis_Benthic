### 1. Map for Study Sites ####
library(sf); library(ggplot2); library(maps); library(ggspatial)
library(dplyr); library(tidyr); library(stringr); library(tidyverse)
library(conflicted)
conflicts_prefer(dplyr::filter)
conflicts_prefer(dplyr::summarise)

rm (list = ls())
getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file")


# North
#FL 25.04, 121.94
#BT 25.13, 121.91
#KI 25.19, 121.78

# Green Island
#GW 22.64, 121.48
#DBS 22.63, 121.49
#GG 22.68, 121.49

# Xiaoliuqiu
#LC 22.34, 120.39
#VN 22.33, 120.35
#DF 22.32, 120.37



## Taiwan ####
# Load the shapefile
taiwan_shapefile <- st_read("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file/gadm41_TWN_0.shp")  # Adjust the path and file name accordingly

# Create a data frame with longitude and latitude coordinates
points_df <- data.frame(
  Longitude = c(121.948, 121.91 , 121.77,
                121.48, 121.49, 121.496,
                120.392, 120.357, 120.373),  # longitudes
  Latitude = c(25.038, 25.14, 25.19,
               22.64, 22.635, 22.68,
               22.344, 22.337, 22.327),  # latitudes
  site = c("FL", "BT", "KI",
           "GW", "DBS", "GG",
           "LC", "VN", "DF"),
  region = c(rep("North", 3), rep("Green Island", 3), rep("Xiaoliuqiu", 3))
)

# Convert to sf object
points_sf <- st_as_sf(points_df, coords = c("Longitude", "Latitude"), crs = 4326)

# Plot the base map with points and labels
base_map <- ggplot(data = taiwan_shapefile) +
  geom_sf(fill = "lightgray", color = "black") +
  theme_minimal() +
  geom_rect(aes(xmin = 121.68, xmax = 122.08, ymin = 24.95, ymax = 25.25), color = "#98c1d9", fill = NA, linetype = "solid") +  # North
  geom_rect(aes(xmin = 121.35, xmax = 121.63, ymin = 22.53, ymax = 22.80), color = "#ffc857", fill = NA, linetype = "solid") +  # Green Island
  geom_rect(aes(xmin = 120.28, xmax = 120.45, ymin = 22.25, ymax = 22.40), color = "#e39695", fill = NA, linetype = "solid") +  # Xiaoliuqiu
  annotation_scale(location = "bl", width_hint = 0.2, height = unit(0.1, "cm"), text_cex = 0.6, line_width = 0.6,
                   pad_x = unit(0.15, "in"), pad_y = unit(0.05, "in")) +  # Custom scale bar
  annotation_north_arrow(location = "tl", which_north = "true", pad_x = unit(0.15, "in"), pad_y = unit(0.05, "in"),
                         style = north_arrow_fancy_orienteering(), height = unit(0.8, "cm"), width = unit(0.8, "cm")) +  # Custom north arrow
  labs(x = "Longitude", y = "Latitude") +
  theme(axis.text.x = element_text(size = 7),  # Adjust x axis text size
        axis.text.y = element_text(size = 7),
        axis.title.x = element_text(size = 9),  # Adjust x axis text size
        axis.title.y = element_text(size = 9))  # Adjust y axis text size


print(base_map)



## North ####
# Create a data frame for North
north_df <- points_df %>% filter(region == "North")
north_df <- north_df %>%
  mutate(vjust = c(-0.9, -1.5, -1.5),
         hjust = c(-0.5, 0.3, 0.3))

# Convert to sf object
north_sf <- st_as_sf(north_df, coords = c("Longitude", "Latitude"), crs = 4326)

# Plot the inset map for North
inset_north <- ggplot(data = north_sf) +
  geom_sf(data = taiwan_shapefile, fill = "lightgray", color = "black") +
  geom_sf(fill = "blue", color = "blue", size = 4) +
  coord_sf(xlim = c(121.68, 122.08), ylim = c(24.95, 25.25)) +
  geom_text(data = north_df, aes(x = Longitude, y = Latitude, label = site, vjust = vjust, hjust = hjust), color = "blue", fontface = "bold", size = 6) +
  annotation_scale(location = "bl", width_hint = 0.3, height = unit(0.3, "cm"), text_cex = 1.0, line_width = 0.6,
                   pad_x = unit(0.15, "in"), pad_y = unit(0.15, "in")) +  # Custom scale bar
  theme_void()

print(inset_north)



## Green Island ####
# Create a data frame for Green Island
green_island_df <- points_df %>% filter(region == "Green Island")
green_island_df <- green_island_df %>%
  mutate(vjust = c(1.6, 1.9, -0.9),
         hjust = c(1.3, 0.5, 0.45))

# Convert to sf object
green_island_sf <- st_as_sf(green_island_df, coords = c("Longitude", "Latitude"), crs = 4326)

# Plot the inset map for Green Island
inset_green_island <- ggplot(data = green_island_sf) +
  geom_sf(data = taiwan_shapefile, fill = "lightgray", color = "black") +
  geom_sf(fill = "blue", color = "blue", size = 4) +
  coord_sf(xlim = c(121.45, 121.52), ylim = c(22.62, 22.69)) +
  geom_text(data = green_island_df, aes(x = Longitude, y = Latitude, label = site, vjust = vjust, hjust = hjust), color = "blue", fontface = "bold", size = 6) +
  annotation_scale(location = "bl", width_hint = 0.3, height = unit(0.3, "cm"), text_cex = 1.0, line_width = 0.6,
                   pad_x = unit(0.05, "in"), pad_y = unit(0.15, "in")) +  # Custom scale bar
  theme_void()

print(inset_green_island)



## Xiaoliuqiu ####
# Create a data frame for Xiaoliuqiu
xiaoliuqiu_df <- points_df %>% filter(region == "Xiaoliuqiu")
xiaoliuqiu_df <- xiaoliuqiu_df %>%
  mutate(vjust = c(-0.9, -1.2, 2.0),
         hjust = c(-0.5, 1.2, -0.3))

# Convert to sf object
xiaoliuqiu_sf <- st_as_sf(xiaoliuqiu_df, coords = c("Longitude", "Latitude"), crs = 4326)

# Plot the inset map for Xiaoliuqiu
inset_xiaoliuqiu <- ggplot(data = xiaoliuqiu_sf) +
  geom_sf(data = taiwan_shapefile, fill = "lightgray", color = "black") +
  geom_sf(fill = "blue", color = "blue", size = 4) +
  coord_sf(xlim = c(120.34, 120.41), ylim = c(22.30, 22.38)) +
  geom_text(data = xiaoliuqiu_df, aes(x = Longitude, y = Latitude, label = site, vjust = vjust, hjust = hjust), color = "blue", fontface = "bold", size = 6) +
  annotation_scale(location = "bl", width_hint = 0.3, height = unit(0.3, "cm"), text_cex = 1.0, line_width = 0.6,
                   pad_x = unit(0.15, "in"), pad_y = unit(0.25, "in")) +  # Custom scale bar
  theme_void()

print(inset_xiaoliuqiu)



######################################################################
### 2. Modify Dataset ####
library(dplyr); library(tidyr); library(stringr); library(tidyverse); library(conflicted)
conflict_prefer("select", 'dplyr')
conflict_prefer("filter", "dplyr")

rm (list = ls())
getwd()

BC.data <- read.csv("benthic_annotations.csv", header=T, sep=",") %>%
  select(Name, Region, Site, Transect, Depth, Label) %>%
  mutate(Depth = ifelse(Depth == 25, 30, Depth), 
         Location = paste0(Site, '_', Depth))


# Labels ####
# 'Genus' represent the most appropriate "Taxon" (most of them are genus level)
BC.labels <- read.csv("labels_OTU.MC.MFG_LN.csv", header=T, sep=",")
BC.labels[c('Genus', 'Morphology')] <- str_split_fixed(BC.labels$OTUs, '_', 2)
BC.labels <- BC.labels[c('Label', 'OTUs', 'Genus', 'Morphology', 'MajorCategory', 'MorphoFunctionalGroup')]
BC.labels <- BC.labels %>%
  mutate(Genus = recode(Genus, 'hard' = 'hard.coral', 
                        'others' = 'others.cca',  'other' = 'other.life',
                        'Red' = 'red.algae',  'Brown' = 'brown.algae',
                        'fishing' = 'fishing.lines',  'fish' = 'fish.net', .default = Genus)) %>%
  mutate(Genus = case_when(
    Morphology == 'spp2_encrusting' ~ 'Merulinidae.spp2',
    Morphology == 'spp2_massive' ~ 'Merulinidae.spp2',
    Morphology == 'spp1_encrusting' ~ 'Merulinidae.spp1',
    Morphology == 'spp1_massive' ~ 'Merulinidae.spp1',
    TRUE ~ Genus)) %>%
  mutate(MajorCategory = recode(MajorCategory, 
                                'Turf' = 'Turf Algae',
                                'CCA' = 'Crustose Coralline Algae',
                                'Black coral' = 'Black Corals',
                                'Gorgonian' = 'Gorgonian Corals',
                                'Soft coral' = 'Soft Corals',
                                'Hard coral' = 'Stony Corals',
                                'Actiniaria' = 'Other Sessile Invertebrates',
                                'Ascidian' = 'Other Sessile Invertebrates',
                                'Hydrozoa' = 'Other Sessile Invertebrates',
                                'Zoanthid' = 'Other Sessile Invertebrates',
                                'Corallimorpharia' = 'Other Sessile Invertebrates',
                                'Fish' = 'Other Life',
                                'Other sessile invertebrates' = 'Other Sessile Invertebrates',
                                'Other mobile invertebrates' = 'Other Mobile Invertebrates',
                                'Stable substrate' = 'Stable Substrate',
                                'Unstable substrate' = 'Unstable Substrate',
                                'Sponge' = 'Sponges',.default = MajorCategory),
         MajorCategory = ifelse(Genus == 'shadow', 'Shadow',
                                ifelse(Genus == 'unknown', 'Unknown', MajorCategory)),
         MajorCategory = recode(MajorCategory, 
                                'Artificial Debris' = 'Other',
                                'Unknown' = 'Other',
                                'TWS' = 'Other',
                                'Shadow' = 'Other',.default = MajorCategory)) %>%
  rename(OTMUs = OTUs,
         OTUs = Genus)


# Full Data (Substrate and Organism, remove Other) ####
BC.fulldata <- BC.data %>% left_join(., BC.labels, by = 'Label') %>%
  mutate(
    #Extract location information from Image name
    Site.code = str_sub(Name, 7,8),
    Site.code = recode(Site.code, 'DB' = 'DBS', .default = Site.code),
    Quadrat = str_sub(Name, -6, -5),
  )
# rearrange the columns' order
BC.fulldata <- BC.fulldata[c('Name', 'Region', 'Site', 'Site.code', 'Depth', 'Location', 'Transect', 'Quadrat', 'Label', 'OTMUs', 'OTUs', 'MajorCategory')]
BC.fulldata <- BC.fulldata %>%
  filter(MajorCategory !='Other')


# Benthic Organism only (remove Substrate) ####
# data only contain the living things
BC.only <- BC.fulldata %>%
  filter(MajorCategory !='Stable Substrate' & MajorCategory !='Unstable Substrate')


# Corals only ####
# data only contain Corals
C.only <- BC.fulldata %>%
  filter(MajorCategory %in% c('Black Corals', 'Gorgonian Corals', 'Soft Corals', 'Stony Corals'))


# Cover ####
BC.cover <- read.csv("benthic_percent_covers.csv", header=T, sep=",")
BC.cover <- head(BC.cover, -1) %>%
  select(-Image.ID, -Annotation.status, -Points) %>%
  pivot_longer(cols = -Image.name, names_to = 'Label', values_to = 'Cover') %>%
  mutate(
    #Extract location information from Image name
    Site = str_sub(Image.name, 7,8),
    Quadrat = str_sub(Image.name, -6, -5),
    Transect = str_sub(Image.name, -9, -8),
    Depth = str_sub(Image.name, -11, -10),
    Region = ifelse(Site == 'KI' | Site == 'BT' | Site == 'FL', 'North', 
                    ifelse(Site == 'LC' | Site == 'VN' | Site == 'DF', 'Xiaoliuqiu', 'Green Island')),
    Depth = ifelse(Depth == '15', '15',
                   ifelse(Depth == '30', '30',
                          ifelse(Depth == '25', '30', '5'))),
    Site = recode(Site, 'DB' = 'DBS', .default = Site))

BC.cover <- BC.cover %>% left_join(., BC.labels, by = 'Label')
BC.cover <- BC.cover %>% select(-Image.name, -OTMUs, -Morphology, -MorphoFunctionalGroup, -Label)
BC.cover <- BC.cover[c('Region', 'Site', 'Depth', 'Transect', 'Quadrat', 'MajorCategory', 'OTUs', 'Cover')]



# BC 3 levels (OTMU -> OTU -> MC) ####
BC.3levels <- BC.fulldata[c('Label', 'OTMUs', 'OTUs', 'MajorCategory')]
BC.note <- read.csv("labels_forNote_LN.csv", header=T, sep=",")
BC.note <- BC.note[c('Label', 'Note')]

BC.3levels <- BC.3levels %>%
  rename(OTU = OTUs,
         OTMU = OTMUs) %>%
  group_by(Label) %>%
  summarise(
    OTMU = names(sort(table(OTMU), decreasing = TRUE)[1]),
    OTU = names(sort(table(OTU), decreasing = TRUE)[1]),
    MajorCategory = names(sort(table(MajorCategory), decreasing = TRUE)[1])
  ) %>%
  ungroup() %>%
  arrange(Label) %>%
  left_join(., BC.note, by = 'Label')


# output the new create csv
# write.csv(BC.labels, 'BC.labels.csv',row.names = F)
# write.csv(BC.fulldata, 'BC.fulldata.csv',row.names = F)
# write.csv(BC.only, 'BC.only.csv',row.names = F)
# write.csv(C.only, 'C.only.csv',row.names = F)
# write.csv(BC.cover, 'BC.cover.csv',row.names = F)
# write.csv(BC.3levels, 'BC.3levels.csv',row.names = T)

# show what Major Category we have in the column: unique(DATA$MajorCategory)



######################################################################
### 3. Alpha Diversity (OTMU richness)---in OTMU level => per Transect OTMUs ####
library(ggplot2); library(dplyr); library(rstatix); library(ggpubr)

rm (list = ls())
getwd()

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



######################################################################
### 4. Bayesian Model ####
# Data manipulation and aggregation / Visualization
library(plyr); library(dplyr); library(ggplot2)
# Bayesian modeling
library(brms); library(rstan); library(stam); library(parallel)
# Model performance evaluation
library(performance)
conflicts_prefer(dplyr::mutate)

rm (list = ls())
getwd()

### Set up the working data ####
BC.cover <- read.csv("BC.cover.csv", header=T, sep=",")


### 4-a. Bayesian Model -- All Coral Cover with Light Attenuation ####
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
# Binomial_Corals.CoverLightA_model <- brm(Coral_points | trials(Tot_Points) ~  Light.Attenuation + (1 | Site), data = Coral_cover2, family = binomial(), control = list(adapt_delta = 0.9, max_treedepth = 11), iter = 4000, warmup = 1000, chains = 2, cores = 2) 
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



### 4-b. Bayesian Model -- Stony Coral Cover with Light Attenuation ####
# to contain only "Stony Coral" we defined
HC.cover <- subset(BC.cover, MajorCategory %in% c("Stony Corals"))
HC.cover <- aggregate (Cover ~ Region + Site + Depth + Transect + Quadrat + MajorCategory + OTUs, HC.cover , sum)

# Modify Coral Cover from per picture 100% to per location (Site + Depth) 100%
HC.cover <- ddply(HC.cover, ~ Region + Site + Depth + MajorCategory + OTUs, function(x){c(Cover = sum(x$Cover))})
HC.cover <- HC.cover %>%
  mutate(Pictures = case_when(
    Depth == 5 ~ 105,  # 105 pictures in 5 m in each site
    Depth == 15 ~ 105, # 105 pictures in 15 m in each site
    Depth == 30 ~ 63)) # 63 pictures in 30 m in each site
HC.cover <- HC.cover %>%
  group_by(Region, Site, Depth, MajorCategory, OTUs) %>%
  reframe(Cover = Cover/Pictures)

HC_data <- HC.cover


### Coral cover profile ####
# Plot considering the effect of Region and Site to see the effect of depth
# Calculate sum for coral cover in function of site and depth for a quick interpretation and plot
HC_cover <- aggregate (Cover ~ Region + Site + Depth, HC_data , sum)

# Transform depth as a Qualitative variable  
HC_cover$Depth <- as.factor(as.character(HC_cover$Depth))
HC_cover$Depth = factor(HC_cover$Depth,levels = c ("5", "15", "30"))

# ggplot with Locations
ggplot(HC_cover, aes(x=Depth, y=Cover)) + 
  geom_boxplot() + geom_point(aes (),size = 1) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=4) + 
  theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))
## mid-domain effect???



### Modify the data for Bayesian model with PAR ####
HC_cover2 <- HC_cover
HC_cover2$Depth <- as.numeric (as.character(HC_cover2$Depth))
HC_cover2$Site <- as.character(HC_cover2$Site)
HC_cover2$Location <- paste(HC_cover2$Site, "_", HC_cover2$Depth)
HC_cover2$RegDep <- paste(HC_cover2$Region, "_", HC_cover2$Depth)
# Add PAR value (use PAR to replace Depth)
PAR.Location <- read.csv("PAR.Location.csv", header=T, sep=",")
HC_cover2 <- HC_cover2 %>%
  left_join(dplyr::select(PAR.Location, Region, Site, Depth, PAR.D.Mean, PAR.S.Mean, Light.Attenuation), by = c("Region", "Site", "Depth"))

# Use a beta distribution to avoid problem of normality
HC_cover2$Cover_Beta <- HC_cover2$Cover/100

HC_cover2$Tot_Points <- 25
HC_cover2$Coral_points <- (HC_cover2$Cover * HC_cover2$Tot_Points) / 100
HC_cover2$Coral_points  <- round(HC_cover2$Coral_points,0) #calculate how many point fall on coral to the nearest integer
HC_cover2$NonCoral_points <- HC_cover2$Tot_Points - HC_cover2$Coral_points 

HC_cover2$Proportion = HC_cover2$Coral_points / (HC_cover2$Coral_points + HC_cover2$NonCoral_points)
HC_cover2$Proportion [HC_cover2$Proportion == 0] <-  .001 # Otherwise beta distribution with 0 it does not work



# Bayesian modelling ####
# Binomial_HC.CoverLightA_model <- brm(Coral_points | trials(Tot_Points) ~  Light.Attenuation + (1 | Site), data = HC_cover2, family = binomial(), control = list(adapt_delta = 0.9, max_treedepth = 11), iter = 4000, warmup = 1000, chains = 2, cores = 2) 
# save(Binomial_HC.CoverLightA_model, file="/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Binomial_HC.CoverLightA_model_Diversity.RData")
load("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file/Binomial_HC.CoverLightA_model_Diversity.RData")


# Summary & Performance evaluation ####
summary (Binomial_HC.CoverLightA_model)
bayes_R2(Binomial_HC.CoverLightA_model)
plot(Binomial_HC.CoverLightA_model)

# "pp_check" help to evaluate how well the Bayesian model fits the data by comparing the observed data with data simulated from the posterior distribution of the model parameters
pp_check(Binomial_HC.CoverLightA_model, type = "scatter_avg") # Not structured data
bayes_R2(Binomial_HC.CoverLightA_model) # 
r2_bayes(Binomial_HC.CoverLightA_model)

HC.me_null <- conditional_effects(Binomial_HC.CoverLightA_model, nsamples = 1000, probs = c(0.05, 0.95), spaghetti = F) # Default is 0.95
plot(HC.me_null, ask = FALSE, points = F) # Probability scale!


# Plot Bayesian Model ####
# Extract the data from HC.me_null
HC_effects_df <- as.data.frame(HC.me_null$Light.Attenuation)
# Predictive plot with confidence interval
bayes_r2_val <- bayes_R2(Binomial_HC.CoverLightA_model)[1]

# Point Color by Region
HC_cover2$Region <- factor(HC_cover2$Region, levels = c("North", "Green Island", "Xiaoliuqiu"))

plot2 <- ggplot(HC_effects_df, aes(x = Light.Attenuation, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "lightgray", alpha = 0.4) +  # Confidence interval
  geom_line(linewidth = 1) +  # Main line
  geom_point(data = HC_cover2, aes(x = Light.Attenuation, y = Cover_Beta, color = Region), size = 2) +  # Real data points
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
    axis.title.x = element_blank(),      # Removes x-axis title
    axis.text.x = element_blank(),        # Removes x-axis text/ticks
    legend.position = "none"  # Removes the legend
  ) +
  annotate("text", x = 0.3, y = 0.85, label = paste("Bayesian R² : ", round(bayes_r2_val, 3)), size = 3) +
  annotate("text", x = 0.3, y = 0.95, label = paste("Stony Corals"), size = 4, fontface = "bold")

plot2



### 4-c. Bayesian Model -- Soft Coral Cover with Light Attenuation ####
# to contain only "Soft Coral" we defined
SC.cover <- subset(BC.cover, MajorCategory %in% c("Soft Corals"))
SC.cover <- aggregate (Cover ~ Region + Site + Depth + Transect + Quadrat + MajorCategory + OTUs, SC.cover , sum)

# Modify Coral Cover from per picture 100% to per location (Site + Depth) 100%
SC.cover <- ddply(SC.cover, ~ Region + Site + Depth + MajorCategory + OTUs, function(x){c(Cover = sum(x$Cover))})
SC.cover <- SC.cover %>%
  mutate(Pictures = case_when(
    Depth == 5 ~ 105,  # 105 pictures in 5 m in each site
    Depth == 15 ~ 105, # 105 pictures in 15 m in each site
    Depth == 30 ~ 63)) # 63 pictures in 30 m in each site
SC.cover <- SC.cover %>%
  group_by(Region, Site, Depth, MajorCategory, OTUs) %>%
  reframe(Cover = Cover/Pictures)

SC_data <- SC.cover


### Coral cover profile ####
# Plot considering the effect of Region and Site to see the effect of depth
# Calculate sum for coral cover in function of site and depth for a quick interpretation and plot
SC_cover <- aggregate (Cover ~ Region + Site + Depth, SC_data , sum)

# Transform depth as a Qualitative variable  
SC_cover$Depth <- as.factor(as.character(SC_cover$Depth))
SC_cover$Depth = factor(SC_cover$Depth,levels = c ("5", "15", "30"))

# ggplot with Locations
ggplot(SC_cover, aes(x=Depth, y=Cover)) + 
  geom_boxplot() + geom_point(aes (),size = 1) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=4) + 
  theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))
## mid-domain effect???



### Modify the data for Bayesian model with PAR ####
SC_cover2 <- SC_cover
SC_cover2$Depth <- as.numeric (as.character(SC_cover2$Depth))
SC_cover2$Site <- as.character(SC_cover2$Site)
SC_cover2$Location <- paste(SC_cover2$Site, "_", SC_cover2$Depth)
SC_cover2$RegDep <- paste(SC_cover2$Region, "_", SC_cover2$Depth)
# Add PAR value (use PAR to replace Depth)
PAR.Location <- read.csv("PAR.Location.csv", header=T, sep=",")
SC_cover2 <- SC_cover2 %>%
  left_join(dplyr::select(PAR.Location, Region, Site, Depth, PAR.D.Mean, PAR.S.Mean, Light.Attenuation), by = c("Region", "Site", "Depth"))

# Use a beta distribution to avoid problem of normality
SC_cover2$Cover_Beta <- SC_cover2$Cover/100

SC_cover2$Tot_Points <- 25
SC_cover2$Coral_points <- (SC_cover2$Cover * SC_cover2$Tot_Points) / 100
SC_cover2$Coral_points  <- round(SC_cover2$Coral_points,0) #calculate how many point fall on coral to the nearest integer
SC_cover2$NonCoral_points <- SC_cover2$Tot_Points - SC_cover2$Coral_points 

SC_cover2$Proportion = SC_cover2$Coral_points / (SC_cover2$Coral_points + SC_cover2$NonCoral_points)
SC_cover2$Proportion [SC_cover2$Proportion == 0] <-  .001 # Otherwise beta distribution with 0 it does not work



# Bayesian modelling ####
# Binomial_SC.CoverLightA_model <- brm(Coral_points | trials(Tot_Points) ~  Light.Attenuation + (1 | Site), data = SC_cover2, family = binomial(), control = list(adapt_delta = 0.9, max_treedepth = 11), iter = 4000, warmup = 1000, chains = 2, cores = 2) 
# save(Binomial_SC.CoverLightA_model, file="/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Binomial_SC.CoverLightA_model_Diversity.RData")
load("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file/Binomial_SC.CoverLightA_model_Diversity.RData")


# Summary & Performance evaluation ####
summary (Binomial_SC.CoverLightA_model)
bayes_R2(Binomial_SC.CoverLightA_model)
plot(Binomial_SC.CoverLightA_model)

# "pp_check" help to evaluate how well the Bayesian model fits the data by comparing the observed data with data simulated from the posterior distribution of the model parameters
pp_check(Binomial_SC.CoverLightA_model, type = "scatter_avg") # Not structured data
bayes_R2(Binomial_SC.CoverLightA_model) # 
r2_bayes(Binomial_SC.CoverLightA_model)

SC.me_null <- conditional_effects(Binomial_SC.CoverLightA_model, nsamples = 1000, probs = c(0.05, 0.95), spaghetti = F) # Default is 0.95
plot(SC.me_null, ask = FALSE, points = F) # Probability scale!


# Plot Bayesian Model ####
# Extract the data from SC.me_null
SC_effects_df <- as.data.frame(SC.me_null$Light.Attenuation)
# Predictive plot with confidence interval
bayes_r2_val <- bayes_R2(Binomial_SC.CoverLightA_model)[1]

# Point Color by Region
SC_cover2$Region <- factor(SC_cover2$Region, levels = c("North", "Green Island", "Xiaoliuqiu"))

plot3 <- ggplot(SC_effects_df, aes(x = Light.Attenuation, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "lightgray", alpha = 0.4) +  # Confidence interval
  geom_line(linewidth = 1) +  # Main line
  geom_point(data = SC_cover2, aes(x = Light.Attenuation, y = Cover_Beta, color = Region), size = 2) +  # Real data points
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
    legend.position = "none"  # Removes the legend
  ) +
  annotate("text", x = 0.3, y = 0.85, label = paste("Bayesian R² : ", round(bayes_r2_val, 3)), size = 3) +
  annotate("text", x = 0.3, y = 0.95, label = paste("Soft Corals"), size = 4, fontface = "bold")

plot3



### 4-d. Bayesian Model -- Black Coral & Gorgonian Coral Cover with Light Attenuation ####
# to contain only "Black Coral & Gorgonian Coral" we defined
BG.cover <- subset(BC.cover, MajorCategory %in% c("Black Corals", "Gorgonian Corals"))
BG.cover <- aggregate (Cover ~ Region + Site + Depth + Transect + Quadrat + MajorCategory + OTUs, BG.cover , sum)

# Modify Coral Cover from per picture 100% to per location (Site + Depth) 100%
BG.cover <- ddply(BG.cover, ~ Region + Site + Depth + MajorCategory + OTUs, function(x){c(Cover = sum(x$Cover))})
BG.cover <- BG.cover %>%
  mutate(Pictures = case_when(
    Depth == 5 ~ 105,  # 105 pictures in 5 m in each site
    Depth == 15 ~ 105, # 105 pictures in 15 m in each site
    Depth == 30 ~ 63)) # 63 pictures in 30 m in each site
BG.cover <- BG.cover %>%
  group_by(Region, Site, Depth, MajorCategory, OTUs) %>%
  reframe(Cover = Cover/Pictures)

BG_data <- BG.cover


### Coral cover profile ####
# Plot considering the effect of Region and Site to see the effect of depth
# Calculate sum for coral cover in function of site and depth for a quick interpretation and plot
BG_cover <- aggregate (Cover ~ Region + Site + Depth, BG_data , sum)

# Transform depth as a Qualitative variable  
BG_cover$Depth <- as.factor(as.character(BG_cover$Depth))
BG_cover$Depth = factor(BG_cover$Depth,levels = c ("5", "15", "30"))

# ggplot with Locations
ggplot(BG_cover, aes(x=Depth, y=Cover)) + 
  geom_boxplot() + geom_point(aes (),size = 1) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=4) + 
  theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))
## mid-domain effect???



### Modify the data for Bayesian model with PAR ####
BG_cover2 <- BG_cover
BG_cover2$Depth <- as.numeric (as.character(BG_cover2$Depth))
BG_cover2$Site <- as.character(BG_cover2$Site)
BG_cover2$Location <- paste(BG_cover2$Site, "_", BG_cover2$Depth)
BG_cover2$RegDep <- paste(BG_cover2$Region, "_", BG_cover2$Depth)
# Add PAR value (use PAR to replace Depth)
PAR.Location <- read.csv("PAR.Location.csv", header=T, sep=",")
BG_cover2 <- BG_cover2 %>%
  left_join(dplyr::select(PAR.Location, Region, Site, Depth, PAR.D.Mean, PAR.S.Mean, Light.Attenuation), by = c("Region", "Site", "Depth"))

# Use a beta distribution to avoid problem of normality
BG_cover2$Cover_Beta <- BG_cover2$Cover/100

BG_cover2$Tot_Points <- 25
BG_cover2$Coral_points <- (BG_cover2$Cover * BG_cover2$Tot_Points) / 100
BG_cover2$Coral_points  <- round(BG_cover2$Coral_points,0) #calculate how many point fall on coral to the nearest integer
BG_cover2$NonCoral_points <- BG_cover2$Tot_Points - BG_cover2$Coral_points 

BG_cover2$Proportion = BG_cover2$Coral_points / (BG_cover2$Coral_points + BG_cover2$NonCoral_points)
BG_cover2$Proportion [BG_cover2$Proportion == 0] <-  .001 # Otherwise beta distribution with 0 it does not work



# Bayesian modelling ####
# Binomial_BG.CoverLightA_model <- brm(Coral_points | trials(Tot_Points) ~  Light.Attenuation + (1 | Site), data = BG_cover2, family = binomial(), control = list(adapt_delta = 0.9, max_treedepth = 11), iter = 4000, warmup = 1000, chains = 2, cores = 2) 
# save(Binomial_BG.CoverLightA_model, file="/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Binomial_BG.CoverLightA_model_Diversity.RData")
load("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file/Binomial_BG.CoverLightA_model_Diversity.RData")


# Summary & Performance evaluation ####
summary (Binomial_BG.CoverLightA_model)
bayes_R2(Binomial_BG.CoverLightA_model)
plot(Binomial_BG.CoverLightA_model)

# "pp_check" help to evaluate how well the Bayesian model fits the data by comparing the observed data with data simulated from the posterior distribution of the model parameters
pp_check(Binomial_BG.CoverLightA_model, type = "scatter_avg") # Not structured data
bayes_R2(Binomial_BG.CoverLightA_model) # 
r2_bayes(Binomial_BG.CoverLightA_model)

BG.me_null <- conditional_effects(Binomial_BG.CoverLightA_model, nsamples = 1000, probs = c(0.05, 0.95), spaghetti = F) # Default is 0.95
plot(BG.me_null, ask = FALSE, points = F) # Probability scale!


# Plot Bayesian Model ####
# Extract the data from BG.me_null
BG_effects_df <- as.data.frame(BG.me_null$Light.Attenuation)
# Predictive plot with confidence interval
bayes_r2_val <- bayes_R2(Binomial_BG.CoverLightA_model)[1]

# Point Color by Region
BG_cover2$Region <- factor(BG_cover2$Region, levels = c("North", "Green Island", "Xiaoliuqiu"))

plot4 <- ggplot(BG_effects_df, aes(x = Light.Attenuation, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "lightgray", alpha = 0.4) +  # Confidence interval
  geom_line(linewidth = 1) +  # Main line
  geom_point(data = BG_cover2, aes(x = Light.Attenuation, y = Cover_Beta, color = Region), size = 2) +  # Real data points
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
    axis.title.x = element_blank(),      # Removes x-axis title
    axis.title.y = element_blank(),      # Removes y-axis title
    axis.text.x = element_blank(),       # Removes x-axis text/ticks
    axis.text.y = element_blank(),       # Removes y-axis text/ticks
    legend.position = "none"  # Removes the legend
  ) +
  annotate("text", x = 0.3, y = 0.85, label = paste("Bayesian R² : ", round(bayes_r2_val, 3)), size = 3) +
  annotate("text", x = 0.25, y = 0.95, label = paste("Black Corals & Gorgonian Corals"), size = 4, fontface = "bold")

plot4



### 4-e. Bayesian Model -- Turf Algae Cover with Light Attenuation ####
# to contain only "Turf Algae" we defined
TA.cover <- subset(BC.cover, MajorCategory %in% c("Turf Algae"))
TA.cover <- aggregate (Cover ~ Region + Site + Depth + Transect + Quadrat + MajorCategory + OTUs, TA.cover , sum)

# Modify Coral Cover from per picture 100% to per location (Site + Depth) 100%
TA.cover <- ddply(TA.cover, ~ Region + Site + Depth + MajorCategory + OTUs, function(x){c(Cover = sum(x$Cover))})
TA.cover <- TA.cover %>%
  mutate(Pictures = case_when(
    Depth == 5 ~ 105,  # 105 pictures in 5 m in each site
    Depth == 15 ~ 105, # 105 pictures in 15 m in each site
    Depth == 30 ~ 63)) # 63 pictures in 30 m in each site
TA.cover <- TA.cover %>%
  group_by(Region, Site, Depth, MajorCategory, OTUs) %>%
  reframe(Cover = Cover/Pictures)

TA_data <- TA.cover


### Coral cover profile ####
# Plot considering the effect of Region and Site to see the effect of depth
# Calculate sum for coral cover in function of site and depth for a quick interpretation and plot
TA_cover <- aggregate (Cover ~ Region + Site + Depth, TA_data , sum)

# Transform depth as a Qualitative variable  
TA_cover$Depth <- as.factor(as.character(TA_cover$Depth))
TA_cover$Depth = factor(TA_cover$Depth,levels = c ("5", "15", "30"))

# ggplot with Locations
ggplot(TA_cover, aes(x=Depth, y=Cover)) + 
  geom_boxplot() + geom_point(aes (),size = 1) + stat_summary(fun=mean, geom="point", shape=18, color="red", size=4) + 
  theme_bw()  + ylab ("Coral cover (%)") + xlab ("Depth (m)") +
  theme(plot.title = element_text(hjust=0.5, size=12, face="bold"),
        axis.text = element_text(size=10, colour="black"),
        axis.title = element_text(size=11, face="bold", colour="black"))
## mid-domain effect???



### Modify the data for Bayesian model with PAR ####
TA_cover2 <- TA_cover
TA_cover2$Depth <- as.numeric (as.character(TA_cover2$Depth))
TA_cover2$Site <- as.character(TA_cover2$Site)
TA_cover2$Location <- paste(TA_cover2$Site, "_", TA_cover2$Depth)
TA_cover2$RegDep <- paste(TA_cover2$Region, "_", TA_cover2$Depth)
# Add PAR value (use PAR to replace Depth)
PAR.Location <- read.csv("PAR.Location.csv", header=T, sep=",")
TA_cover2 <- TA_cover2 %>%
  left_join(dplyr::select(PAR.Location, Region, Site, Depth, PAR.D.Mean, PAR.S.Mean, Light.Attenuation), by = c("Region", "Site", "Depth"))

# Use a beta distribution to avoid problem of normality
TA_cover2$Cover_Beta <- TA_cover2$Cover/100

TA_cover2$Tot_Points <- 25
TA_cover2$Coral_points <- (TA_cover2$Cover * TA_cover2$Tot_Points) / 100
TA_cover2$Coral_points  <- round(TA_cover2$Coral_points,0) #calculate how many point fall on coral to the nearest integer
TA_cover2$NonCoral_points <- TA_cover2$Tot_Points - TA_cover2$Coral_points 

TA_cover2$Proportion = TA_cover2$Coral_points / (TA_cover2$Coral_points + TA_cover2$NonCoral_points)
TA_cover2$Proportion [TA_cover2$Proportion == 0] <-  .001 # Otherwise beta distribution with 0 it does not work



# Bayesian modelling ####
# Binomial_TA.CoverLightA_model <- brm(Coral_points | trials(Tot_Points) ~  Light.Attenuation + (1 | Site), data = TA_cover2, family = binomial(), control = list(adapt_delta = 0.9, max_treedepth = 11), iter = 4000, warmup = 1000, chains = 2, cores = 2) 
# save(Binomial_TA.CoverLightA_model, file="/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Binomial_TA.CoverLightA_model_Diversity.RData")
load("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file/Binomial_TA.CoverLightA_model_Diversity.RData")


# Summary & Performance evaluation ####
summary (Binomial_TA.CoverLightA_model)
bayes_R2(Binomial_TA.CoverLightA_model)
plot(Binomial_TA.CoverLightA_model)

# "pp_check" help to evaluate how well the Bayesian model fits the data by comparing the observed data with data simulated from the posterior distribution of the model parameters
pp_check(Binomial_TA.CoverLightA_model, type = "scatter_avg") # Not structured data
bayes_R2(Binomial_TA.CoverLightA_model) # 
r2_bayes(Binomial_TA.CoverLightA_model)

TA.me_null <- conditional_effects(Binomial_TA.CoverLightA_model, nsamples = 1000, probs = c(0.05, 0.95), spaghetti = F) # Default is 0.95
plot(TA.me_null, ask = FALSE, points = F) # Probability scale!


# Plot Bayesian Model ####
# Extract the data from TA.me_null
TA_effects_df <- as.data.frame(TA.me_null$Light.Attenuation)
# Predictive plot with confidence interval
bayes_r2_val <- bayes_R2(Binomial_TA.CoverLightA_model)[1]

# Point Color by Region
TA_cover2$Region <- factor(TA_cover2$Region, levels = c("North", "Green Island", "Xiaoliuqiu"))

plot5 <- ggplot(TA_effects_df, aes(x = Light.Attenuation, y = estimate__)) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "lightgray", alpha = 0.4) +  # Confidence interval
  geom_line(linewidth = 1) +  # Main line
  geom_point(data = TA_cover2, aes(x = Light.Attenuation, y = Cover_Beta, color = Region), size = 2) +  # Real data points
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
    axis.title.x = element_blank(),      # Removes x-axis title
    axis.title.y = element_blank(),      # Removes y-axis title
    axis.text.x = element_blank(),       # Removes x-axis text/ticks
    axis.text.y = element_blank(),       # Removes y-axis text/ticks
    legend.position = "none"  # Removes the legend
  ) +
  annotate("text", x = 0.3, y = 0.85, label = paste("Bayesian R² : ", round(bayes_r2_val, 3)), size = 3) +
  annotate("text", x = 0.3, y = 0.95, label = paste("Turf Algae"), size = 4, fontface = "bold")

plot5



### 4-f. Bayesian Model -- CCA Cover with Light Attenuation ####
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
# Binomial_CCA.CoverLightA_model <- brm(Coral_points | trials(Tot_Points) ~  Light.Attenuation + (1 | Site), data = CCA_cover2, family = binomial(), control = list(adapt_delta = 0.9, max_treedepth = 11), iter = 4000, warmup = 1000, chains = 2, cores = 2) 
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
plotA <- plot1 + plot4
plotB <- plot2 + plot5
plotC <- plot3 + plot6
plotA / plotB / plotC
# Combine all plots and add legend at the bottom
combined_plot <- (plotA / plotB / plotC) + plot_layout(guides = "collect")
combined_plot



######################################################################
### 5. Beta Diversity: Turnover & Nestedness of BC & CA ####
## Beta Diversity: Turnover & Nestedness => heatmap---in OTMU level
library(dplyr); library(tidyr); library(tibble)
library(betapart); library(graphics)
conflicts_prefer(dplyr::summarize)
conflicts_prefer(dplyr::filter)

rm (list = ls())
getwd()

### *all Benthic Community data---in OTMU level ####
BC.only <- read.csv("BC.only.csv", header = T, sep = ",")

OTMU_bc.sum <- BC.only %>%
  group_by(Region, Location, Depth, MajorCategory, OTMUs) %>%
  summarise(Number = length(OTMUs))
OTMU_bc.sum$Region <- factor(OTMU_bc.sum$Region, 
                             levels = c('North', 'Green Island','Xiaoliuqiu'))
OTMU_bc.sum$Depth <- factor(OTMU_bc.sum$Depth, 
                            levels = c('5', '15','30'))
OTMU_bc.sum <- OTMU_bc.sum %>%
  mutate(Reg.Dep = factor(paste(Region, Depth, sep = "_"),
                          levels = c("North_5", "North_15", "North_30",
                                     "Green Island_5", "Green Island_15", "Green Island_30",
                                     "Xiaoliuqiu_5", "Xiaoliuqiu_15", "Xiaoliuqiu_30")))

beta.BC <- OTMU_bc.sum %>%
  group_by(Reg.Dep, OTMUs) %>%
  summarise(Abundance = sum(Number), .groups = 'drop') %>% # remove grouping
  pivot_wider(names_from = OTMUs, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  column_to_rownames(var = "Reg.Dep")

# Convert to numeric matrix
beta.BC <- as.matrix(beta.BC)
beta.BC <- beta.BC[levels(OTMU_bc.sum$Reg.Dep), ]
# Convert to presence/absence matrix
beta_jar_mat <- ifelse(beta.BC > 1, 1, 0)

# Calculate beta diversity
beta.multi_res <- beta.multi(beta_jar_mat, index.family = "jaccard")
beta_location <- beta.pair(beta_jar_mat, index.family = "jaccard")

# Define a color palette: adjust the colors and the number of levels
my_palette <- colorRampPalette(c("white", "yellow", "orange", "red"))(200)
rev_my_palette <- rev(my_palette)


## 5-a. Beta Diversity of BC ####
beta_location_ac <- beta_location$beta.jac
dst <- data.matrix(beta_location_ac)

dst_ac_nor <- dst[1:9,1:9] # order => 1:3->N / 4:6->GI / 7:9->XLQ
x_labels <- c("NT 5", "NT 15", "NT 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("NT 5", "NT 15", "NT 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")

plot_BC_beta <- function() {
  # Full plot of the heatmap
  image(1:dim(dst_ac_nor)[1], 1:dim(dst_ac_nor)[2], dst_ac_nor, col = my_palette, axes = F, xlab = "", ylab = "")
  axis(1, at = 1:dim(dst_ac_nor)[1], labels = x_labels, cex.axis = 1.0, font.axis = 2, las = 1)
  axis(2, at = 1:dim(dst_ac_nor)[2], labels = y_labels, cex.axis = 1.0, font.axis = 2, las = 1)
  text(expand.grid(1:dim(dst_ac_nor)[1], 1:dim(dst_ac_nor)[2]), sprintf("%0.2f", dst_ac_nor), cex = 1.5, font = 2)
  
  # Extract the upper triangular part
  dst_ac_nor_upper <- dst_ac_nor
  dst_ac_nor_upper[lower.tri(dst_ac_nor_upper)] <- NA
  
  # Create the heatmap plot for the upper triangular part
  image(1:dim(dst_ac_nor)[1], 1:dim(dst_ac_nor)[2], dst_ac_nor_upper, col = my_palette, axes = FALSE, xlab = "", ylab = "")
  axis(3, at = 1:dim(dst_ac_nor)[1], labels = x_labels, cex.axis = 1.0, font.axis = 2, las = 1)
  axis(2, at = 1:dim(dst_ac_nor)[2], labels = y_labels, cex.axis = 1.0, font.axis = 2, las = 1)
  
  # Add text annotations for the upper triangular part
  text_coords <- expand.grid(1:dim(dst_ac_nor)[1], 1:dim(dst_ac_nor)[2])
  text_coords <- text_coords[upper.tri(dst_ac_nor), ]
  text(text_coords, labels = sprintf("%0.2f", dst_ac_nor_upper[upper.tri(dst_ac_nor_upper)]), cex = 1.5, font = 2)
}
plot_BC_beta()
graphics.off()
  
  
## Turnover of BC ####
beta_location_tu <- beta_location$beta.jtu
dst <- data.matrix(beta_location_tu)

dst_tu_nor <- dst[1:9, 1:9]  # order => 1:3->N / 4:6->GI / 7:9->XLQ
x_labels <- c("NT 5", "NT 15", "NT 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("NT 5", "NT 15", "NT 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")


## Nestedness of BC ####
beta_location_ne <- beta_location$beta.jne
dst <- data.matrix(beta_location_ne)

dst_ne_nor <- dst[1:9,1:9] # order => 1:3->N / 4:6->GI / 7:9->XLQ
x_labels <- c("NT 5", "NT 15", "NT 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("NT 5", "NT 15", "NT 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")


## *conbine Turnover & Nestedness of BC map together ####
# Extract the upper triangular part for Turnover
dst_tu_nor_upper <- dst_tu_nor
dst_tu_nor_upper[lower.tri(dst_tu_nor_upper)] <- NA
# Extract the lower triangular part for Nestedness
dst_ne_nor_lower <- dst_ne_nor
dst_ne_nor_lower[upper.tri(dst_ne_nor_lower)] <- NA

# Combine the upper triangular and lower triangular parts
combined_matrix <- dst_tu_nor_upper
combined_matrix[lower.tri(combined_matrix)] <- dst_ne_nor_lower[lower.tri(dst_ne_nor_lower)]

# Set layout for the heatmap and the legend
layout(matrix(c(1, 2), nrow = 1), widths = c(6, 1)) # Adjust widths as needed
# Create the heatmap plot for the combined matrix
par(mar = c(5, 5, 2, 2)) # Adjust margins as needed
plot_BC_TN <- function() {
  image(1:dim(combined_matrix)[1], 1:dim(combined_matrix)[2], combined_matrix, col = my_palette, axes = FALSE, xlab = "", ylab = "")
  axis(1, at = 1:dim(combined_matrix)[1], labels = x_labels, cex.axis = 1.0, font.axis = 2, las = 1)
  axis(2, at = 1:dim(combined_matrix)[2], labels = y_labels, cex.axis = 1.0, font.axis = 2, las = 1)
  
  # Add text annotations for the upper triangular part (Turnover)
  text_coords_upper <- expand.grid(1:dim(combined_matrix)[1], 1:dim(combined_matrix)[2])
  text_coords_upper <- text_coords_upper[upper.tri(combined_matrix), ]
  text(text_coords_upper, labels = sprintf("%0.2f", combined_matrix[upper.tri(combined_matrix)]), cex = 1.5, font = 2)
  # Add text annotations for the lower triangular part (Nestedness)
  text_coords_lower <- expand.grid(1:dim(combined_matrix)[1], 1:dim(combined_matrix)[2])
  text_coords_lower <- text_coords_lower[lower.tri(combined_matrix), ]
  text(text_coords_lower, labels = sprintf("%0.2f", combined_matrix[lower.tri(combined_matrix)]), cex = 1.5, font = 2)
  
  # Add a bold line from top-right to bottom-left
  segments(x0 = dim(combined_matrix)[1] + 0.5, y0 = dim(combined_matrix)[2] + 0.5, x1 = 0.5, y1 = 0.5, lwd = 4, col = "black")
  # dev.off()
}
plot_BC_TN()
graphics.off()
  
  
### *only Coral Assemblage data---in OTMU level ####
Coral.only <- BC.only %>%
  filter(MajorCategory =='Black Corals'| MajorCategory =='Gorgonian Corals'| MajorCategory =='Soft Corals'| MajorCategory =='Stony Corals')

OTMU_coral.sum <- Coral.only %>%
  group_by(Region, Location, Depth, MajorCategory, OTMUs) %>%
  summarise(Number = length(OTMUs))
OTMU_coral.sum$Region <- factor(OTMU_coral.sum$Region, 
                                levels = c('North', 'Green Island','Xiaoliuqiu'))
OTMU_coral.sum$Depth <- factor(OTMU_coral.sum$Depth, 
                               levels = c('5', '15','30'))
OTMU_coral.sum <- OTMU_coral.sum %>%
  mutate(Reg.Dep = factor(paste(Region, Depth, sep = "_"),
                          levels = c("North_5", "North_15", "North_30",
                                     "Green Island_5", "Green Island_15", "Green Island_30",
                                     "Xiaoliuqiu_5", "Xiaoliuqiu_15", "Xiaoliuqiu_30")))

beta.CA <- OTMU_coral.sum %>%
  group_by(Reg.Dep, OTMUs) %>%
  summarise(Abundance = sum(Number), .groups = 'drop') %>% # remove grouping
  pivot_wider(names_from = OTMUs, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  column_to_rownames(var = "Reg.Dep")

# Convert to numeric matrix
beta.CA <- as.matrix(beta.CA)
beta.CA <- beta.CA[levels(OTMU_coral.sum$Reg.Dep), ]
# Convert to presence/absence matrix
beta_jar_mat.ca <- ifelse(beta.CA > 1, 1, 0)

# Calculate beta diversity
beta.multi_res.ca <- beta.multi(beta_jar_mat.ca, index.family = "jaccard")
beta_location.ca <- beta.pair(beta_jar_mat.ca, index.family = "jaccard")

# Define a color palette: adjust the colors and the number of levels
my_palette <- colorRampPalette(c("white", "yellow", "orange", "red"))(200)
rev_my_palette <- rev(my_palette)

## 5-b. Beta Diversity of CA ####
beta_location.ca_ac <- beta_location.ca$beta.jac
dst <- data.matrix(beta_location.ca_ac)

dst_ac_nor <- dst[1:9,1:9] # order => 1:3->N / 4:6->GI / 7:9->XLQ
x_labels <- c("NT 5", "NT 15", "NT 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("NT 5", "NT 15", "NT 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")


plot_CA_beta <- function() {
  # Full plot of the heatmap
  image(1:dim(dst_ac_nor)[1], 1:dim(dst_ac_nor)[2], dst_ac_nor, col = my_palette, axes = F, xlab = "", ylab = "")
  axis(1, at = 1:dim(dst_ac_nor)[1], labels = x_labels, cex.axis = 1.0, font.axis = 2, las = 1)
  axis(2, at = 1:dim(dst_ac_nor)[2], labels = y_labels, cex.axis = 1.0, font.axis = 2, las = 1)
  text(expand.grid(1:dim(dst_ac_nor)[1], 1:dim(dst_ac_nor)[2]), sprintf("%0.2f", dst_ac_nor), cex = 1.5, font = 2)
  
  # Extract the upper triangular part
  dst_ac_nor_upper <- dst_ac_nor
  dst_ac_nor_upper[lower.tri(dst_ac_nor_upper)] <- NA
  
  # Create the heatmap plot for the upper triangular part
  image(1:dim(dst_ac_nor)[1], 1:dim(dst_ac_nor)[2], dst_ac_nor_upper, col = my_palette, axes = FALSE, xlab = "", ylab = "")
  axis(3, at = 1:dim(dst_ac_nor)[1], labels = x_labels, cex.axis = 1.0, font.axis = 2, las = 1)
  axis(2, at = 1:dim(dst_ac_nor)[2], labels = y_labels, cex.axis = 1.0, font.axis = 2, las = 1)
  
  # Add text annotations for the upper triangular part
  text_coords <- expand.grid(1:dim(dst_ac_nor)[1], 1:dim(dst_ac_nor)[2])
  text_coords <- text_coords[upper.tri(dst_ac_nor), ]
  text(text_coords, labels = sprintf("%0.2f", dst_ac_nor_upper[upper.tri(dst_ac_nor_upper)]), cex = 1.5, font = 2)
}
plot_CA_beta()
graphics.off()
  
  
## Turnover of CA ####
beta_location.ca_tu <- beta_location.ca$beta.jtu
dst <- data.matrix(beta_location.ca_tu)

dst_tu_nor <- dst[1:9, 1:9]  # order => 1:3->N / 4:6->GI / 7:9->XLQ
x_labels <- c("NT 5", "NT 15", "NT 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("NT 5", "NT 15", "NT 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")


## Nestedness of CA ####
beta_location.ca_ne <- beta_location.ca$beta.jne
dst <- data.matrix(beta_location.ca_ne)

dst_ne_nor <- dst[1:9,1:9] # order => 1:3->N / 4:6->GI / 7:9->XLQ
x_labels <- c("NT 5", "NT 15", "NT 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("NT 5", "NT 15", "NT 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")


## *conbine Turnover & Nestedness of CA map together ####
# Extract the upper triangular part for Turnover
dst_tu_nor_upper <- dst_tu_nor
dst_tu_nor_upper[lower.tri(dst_tu_nor_upper)] <- NA
# Extract the lower triangular part for Nestedness
dst_ne_nor_lower <- dst_ne_nor
dst_ne_nor_lower[upper.tri(dst_ne_nor_lower)] <- NA

# Combine the upper triangular and lower triangular parts
combined_matrix <- dst_tu_nor_upper
combined_matrix[lower.tri(combined_matrix)] <- dst_ne_nor_lower[lower.tri(dst_ne_nor_lower)]

# Set layout for the heatmap and the legend
layout(matrix(c(1, 2), nrow = 1), widths = c(6, 1)) # Adjust widths as needed
# Create the heatmap plot for the combined matrix
par(mar = c(5, 5, 2, 2)) # Adjust margins as needed
plot_CA_TN <- function() {
  image(1:dim(combined_matrix)[1], 1:dim(combined_matrix)[2], combined_matrix, col = my_palette, axes = FALSE, xlab = "", ylab = "")
  axis(1, at = 1:dim(combined_matrix)[1], labels = x_labels, cex.axis = 1.0, font.axis = 2, las = 1)
  axis(2, at = 1:dim(combined_matrix)[2], labels = y_labels, cex.axis = 1.0, font.axis = 2, las = 1)
  
  # Add text annotations for the upper triangular part (Turnover)
  text_coords_upper <- expand.grid(1:dim(combined_matrix)[1], 1:dim(combined_matrix)[2])
  text_coords_upper <- text_coords_upper[upper.tri(combined_matrix), ]
  text(text_coords_upper, labels = sprintf("%0.2f", combined_matrix[upper.tri(combined_matrix)]), cex = 1.5, font = 2)
  # Add text annotations for the lower triangular part (Nestedness)
  text_coords_lower <- expand.grid(1:dim(combined_matrix)[1], 1:dim(combined_matrix)[2])
  text_coords_lower <- text_coords_lower[lower.tri(combined_matrix), ]
  text(text_coords_lower, labels = sprintf("%0.2f", combined_matrix[lower.tri(combined_matrix)]), cex = 1.5, font = 2)
  
  # Add a bold line from top-right to bottom-left
  segments(x0 = dim(combined_matrix)[1] + 0.5, y0 = dim(combined_matrix)[2] + 0.5, x1 = 0.5, y1 = 0.5, lwd = 4, col = "black")
  # dev.off()
}
plot_CA_TN()
graphics.off()



######################################################################
### 6. NMDS / PERMANOVA of BC & CA ####
## beta quantitative => NMDS + bray curtis + permanova + simper (indicative species)
library(dplyr); library(tidyr); library(tibble)
library(vegan); library(labdsv); library(RVAideMemoire)
# Load the pairwiseAdonis package if not already done
# install.packages("devtools")
# devtools::install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

rm (list = ls())
getwd()

### 6-a. Benthic Community NMDS---in OTMU level ####
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



## 6-b. Coral Assemblage NMDS---in OTMU level ####
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



######################################################################
### 7. Light Attenuation / Benthic Composition / Top 5 Coral OTU ####
library(ggplot2); library(dplyr); library(scales); library(tidyr)
conflicts_prefer(dplyr::arrange)

rm (list = ls())
getwd()

### 7-1. Light Attenuation ####
PAR.diving.data <- read.csv("PAR raw data diving.csv", header=T, sep=",")
PAR.surface.data <- read.csv("PAR raw data surface.csv", header=T, sep=",")

## *Light Attenuation Data ####
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
  mutate(Attenuation = D.Mean / S.Mean)
# each location's Attenuation value
# write.csv(PAR.full, 'Light_Locations Attenuation.csv',row.names = T)
PAR.full <- PAR.full %>%
  group_by(Region, Depth) %>%
  summarise(A.SD = sd(Attenuation, na.rm = F),
            Attenuation = mean(Attenuation),
            D.SD = sd(D.Mean, na.rm = F),
            S.SD = sd(S.Mean, na.rm = F),
            D.Mean = mean(D.Mean),
            S.Mean = mean(S.Mean))
PAR.full <- PAR.full[c('Region', 'Depth', 'D.Mean', 'D.SD', 'S.Mean', 'S.SD', 'Attenuation', 'A.SD')]
PAR.full$Region <- factor(PAR.full$Region, levels = c('North', 'Green Island','Xiaoliuqiu'))
PAR.full$Depth <- factor(PAR.full$Depth, levels = c('30', '15','5'))
# write.csv(PAR.full, 'PAR value & Light Attenuation.csv',row.names = T)


## *PAR Region SD ####
PAR.RS.sd <- PAR.full %>%
  group_by(Region) %>%
  summarise(R.SD = sd(S.Mean, na.rm = F),
            R.S.Mean = mean(S.Mean))

PAR.RS <- PAR.full %>%
  group_by(Region) %>%
  summarise(S.SD = sd(S.Mean, na.rm = F),
            S.Mean = mean(S.Mean),
            D.Mean = S.Mean,
            D.SD = S.SD,
            Attenuation = D.Mean / S.Mean)

## *add 0 m PAR data
PAR.RS <- PAR.RS %>%
  left_join(PAR.RS.sd, by = "Region") %>% # Join R.SD from PAR.RS.sd
  mutate(Depth = 0,
         D.Mean = S.Mean,
         D.SD = R.SD,
         S.SD = R.SD,
         A.SD = sd(Attenuation, na.rm = FALSE)) %>%
  select(Region, Depth, D.Mean, D.SD, S.Mean, S.SD, Attenuation, A.SD)
PAR.RS <- PAR.RS[c('Region', 'Depth', 'D.Mean', 'D.SD', 'S.Mean', 'S.SD', 'Attenuation', 'A.SD')]
# Bind the new rows to the original dataset
PAR.full$Depth <- as.numeric(as.character(PAR.full$Depth))
PAR.full2 <- bind_rows(PAR.full, PAR.RS)
PAR.full2$Depth <- factor(PAR.full2$Depth, levels = c('0', '5', '15', '30'))
PAR.full2$Region <- factor(PAR.full2$Region, levels = c('North', 'Green Island', 'Xiaoliuqiu'))
PAR.full2 <- PAR.full2 %>%
  arrange(factor(Region, levels = c('North', 'Green Island', 'Xiaoliuqiu')),
          factor(Depth, levels = c('0', '5', '15', '30')))
# write.csv(PAR.full2, 'PAR value & Light Attenuation (0-30).csv',row.names = F)

## *Attenuation Plot ####
PAR.full3 <- PAR.full2
PAR.full3 <- PAR.full3[c('Region', 'Depth', 'D.Mean', 'D.SD', 'S.Mean', 'S.SD', 'Attenuation', 'A.SD')]
PAR.full3$Depth <- factor(PAR.full3$Depth, levels = c('30', '15', '5', '0'))
PAR.full3$Region <- factor(PAR.full3$Region, levels = c('North', 'Green Island','Xiaoliuqiu'))

PAR <- PAR.full3 %>%
  group_by(Region) %>%
  mutate(xend = lead(Attenuation),
         yend = lead(Depth),
         linetype = factor(row_number()))
#drop_na() %>%

# Convert factor columns to character
PAR <- PAR %>%
  mutate(yend = ifelse(is.na(yend) == T, 0 , yend),
         xend = ifelse(is.na(xend) == T, 0 , xend)) 

# plot
PAR %>%
  ggplot(aes(x = Attenuation, y = Depth, color = Region)) +
  geom_segment(aes(xend = xend, yend = yend, linetype = linetype), linewidth = 0.9) +
  scale_linetype_manual(values = c("dashed", "solid", "solid", "blank")) +
  
  geom_line(aes(color=Region), linewidth=1.0) +
  geom_point(size=1.2, color = "black") +
  geom_errorbarh(aes(xmin = Attenuation - A.SD, 
                     xmax = Attenuation + A.SD, 
                     color = Region), 
                 height = 0.1,
                 size = 0.5) +
  geom_text(aes(label = sprintf("%.1f", D.Mean),
                color = Region,  # Match color to Region
                vjust = case_when(
                  Depth == 0 & Region == "North" ~ -3.5,
                  Depth == 0 & Region == "Green Island" ~ -2.0,
                  Depth == 0 & Region == "Xiaoliuqiu" ~ -0.5,
                  Depth == 5 & Region == "North" ~ -0.9,
                  Depth == 15 & Region == "North" ~ -1.0,
                  Depth == 30 & Region == "North" ~ -0.9,
                  Depth == 5 & Region == "Green Island" ~ 2.0,
                  Depth == 15 & Region == "Green Island" ~ -0.65,
                  Depth == 30 & Region == "Green Island" ~ 2.0,
                  Depth == 5 & Region == "Xiaoliuqiu" ~ -1.7,
                  Depth == 15 & Region == "Xiaoliuqiu" ~ 1.9,
                  Depth == 30 & Region == "Xiaoliuqiu" ~ 2.0),
                hjust = case_when(
                  Depth == 5 & Region == "North" ~ 1.0,
                  Depth == 15 & Region == "North" ~ 1.0,
                  Depth == 30 & Region == "North" ~ 1.1,
                  Depth == 5 & Region == "Green Island" ~ -0.1,
                  Depth == 15 & Region == "Green Island" ~ -0.55,
                  Depth == 30 & Region == "Green Island" ~ 0.8,
                  Depth == 5 & Region == "Xiaoliuqiu" ~ 0.4,
                  Depth == 15 & Region == "Xiaoliuqiu" ~ -0.1,
                  Depth == 30 & Region == "Xiaoliuqiu" ~ -0.3)),
            size = 2.8) +
  scale_color_manual(values=c('#98c1d9','#ffc857','#adc178'),
                     labels = c('North Taiwan', 'Green Island', 'Xiaoliuqiu')) + # change label name "North" to "North Taiwan"
  scale_x_continuous(name = "% of Surface PAR",
                     limits = c(0, 1.00),
                     breaks = seq(0, 1.00, 0.2),
                     labels = seq(0, 100, 20),  # Change labels to whole numbers
                     position = "top") +
  scale_y_discrete(name = "Depth (m)") +
  theme(legend.position="right", 
        plot.title = element_text(hjust = 0.5, size = 20),
        axis.title = element_text(size = 14, face = "bold"),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        legend.text = element_text(size = 12),
        legend.key.size = unit(0.8, 'cm'),
        panel.background = element_blank(),  # Removes the gray background of the plot panel
        panel.grid.major = element_line(color = "#e9ecef"),  # Adds major grid lines
        panel.grid.minor = element_line(color = "#e9ecef"),
        axis.ticks = element_line(color = "black"),  # Keeps the ticks
        axis.line = element_line(color = "black")
  ) +
  guides(linetype = "none")  # Hide the linetype legend only


### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
# PAR data for Bayesian Model ####
PAR.Location <- PAR.diving.mean %>%
  left_join(select(PAR.surface.mean, Site, Depth, S.Mean), by = c("Site", "Depth")) %>%
  mutate(Attenuation = D.Mean / S.Mean) %>%
  group_by(Region, Site, Depth) %>%
  summarise(PAR.D.Mean = mean(D.Mean),
            PAR.S.Mean = mean(S.Mean),
            Light.Attenuation = mean(Attenuation))

PAR.Location$Region <- factor(PAR.Location$Region, levels = c('North', 'Green Island','Xiaoliuqiu'))
PAR.Location$Depth <- factor(PAR.Location$Depth, levels = c('5', '15', '30'))
# write.csv(PAR.Location, 'PAR.Location.csv',row.names = F)
# write.csv(PAR.full, 'PAR.RegDep.csv',row.names = F)
### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###



### 7-2. Benthic Composition---in Major Category level ####
BC.fulldata <- read.csv("BC.fulldata.csv", header=T, sep=",")

# ***calculate % for BC ####
MC_sum.numb <- BC.fulldata %>% group_by(Region,Location,Depth,MajorCategory) %>% 
  summarise(Number=length(MajorCategory)) %>% ungroup()
## show total Major Category in the column
# MC_list <- unique(MC_sum.numb$MajorCategory)
# print(MC_list)

MC_L.sum.numb <- MC_sum.numb %>%
  mutate(RegDep = paste0(Region, '_', Depth)) %>%
  group_by(Region, RegDep, Location, Depth ,MajorCategory) %>%
  summarise(Number = sum(Number))

# Define all MajorCategory levels
all_categories <- c('Turf Algae', 'Macroalgae', 'Crustose Coralline Algae',
                    'Sponges', 'Other Sessile Invertebrates', 'Other Mobile Invertebrates', 'Other Life',
                    'Unstable Substrate', 'Stable Substrate',
                    'Black Corals', 'Gorgonian Corals', 'Soft Corals', 'Stony Corals')
# Create a complete dataset with all combinations of Region, Depth, Location, and MajorCategory
complete_data <- MC_L.sum.numb %>%
  distinct(Region, Depth, Location) %>% # Get all unique combinations of Region, Depth, and Location
  crossing(MajorCategory = all_categories) # Add all categories to each combination
# Join the complete dataset with the original data
MC_L.sum.numb_filled <- complete_data %>%
  left_join(MC_L.sum.numb, by = c("Region", "Depth", "RegDep", "Location", "MajorCategory")) %>%
  mutate(Number = replace_na(Number, 0)) # Replace NA with 0 for missing combinations


MC_BC.percent2 <- MC_L.sum.numb_filled %>%
  group_by(Location) %>%
  reframe(Total = sum(Number)) %>%
  right_join(MC_L.sum.numb_filled, by = c("Location")) %>%
  mutate(Percent = Number / Total) %>%
  group_by(Region, Location, Depth ,MajorCategory) %>%
  summarise(Percent = Percent,
            Percent2 = sprintf("%.1f%%", Percent * 100), 
            .groups = 'drop')
MC_BC.percent3 <- MC_BC.percent2 %>%
  group_by(Region, Depth, MajorCategory) %>%
  summarise(SD = sd(Percent, na.rm = F),
            Percent = mean(Percent),
            Percent2 = Percent * 100,
            Percent3 = sprintf("%.1f%%", Percent * 100),
            SD3 = sprintf("%.1f%%", SD * 100),
            .groups = 'drop')
MC_BC.percent3 <- MC_BC.percent3[c('Region', 'Depth', 'MajorCategory', 'Percent', 'SD', 'Percent2', 'Percent3', 'SD3')]


MC_BC.percent3$MajorCategory <- factor(MC_BC.percent3$MajorCategory, 
                                       levels = c('Turf Algae', 'Macroalgae', 'Crustose Coralline Algae',
                                                  'Sponges', 'Other Sessile Invertebrates', 'Other Mobile Invertebrates', 'Other Life',
                                                  'Unstable Substrate', 'Stable Substrate',
                                                  'Black Corals', 'Gorgonian Corals', 'Soft Corals', 'Stony Corals'))
MC_BC.percent3$Region <- factor(MC_BC.percent3$Region, levels = c('North', 'Green Island','Xiaoliuqiu'))
MC_BC.percent3$Depth <- factor(MC_BC.percent3$Depth, levels = c('5', '15','30'))
# write.csv(MC_BC.percent3, 'BC_MajorCategory.percent.csv', row.names = F)

# BC plot ####
ggplot(data = MC_BC.percent3, mapping = aes(x = Depth, y = Percent2, fill = MajorCategory)) + 
  geom_bar(position = "fill", stat = 'identity') + 
  labs(fill = "Major Category") + 
  ylab("Cover (%)") + 
  xlab("Depth (m)") + 
  scale_fill_manual(values = c("#dde5b6", "#adc178", "#c9cba3",
                               "#a9d6e5", "#61a5c2", "#2a6f97", "#3d5a80",
                               "#adb5bd", "#495057",
                               "#fee8c8", "#fdd49e", "#fdab74", "#fc8d59")) + #"#ced4da", 
  facet_wrap(~Region, scales = "free_x") +  # Keep x-axis ticks across all facets
  scale_y_continuous(labels = seq(0, 100, 25)) +
  theme(
    text = element_text(size = 12),            # Base text size
    axis.title = element_text(size = 14, face = "bold"),      # Axis titles text size
    axis.text = element_text(size = 12),       # Axis labels text size
    strip.text = element_text(size = 13, face = "bold"),      # Facet labels text size
    legend.text = element_text(size = 10),     # Legend text size
    legend.title = element_text(size = 12),     # Legend title text size
    legend.key.size = unit(0.6, 'cm'),  # Adjusts the overall size of the legend keys
    panel.background = element_blank(),  # Removes the gray background of the plot panel
    plot.background = element_blank(),    # Removes the gray background of the entire plot
    axis.ticks = element_line(color = "black"),  # Keeps the ticks
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Adds black borders around each panel
    strip.background = element_rect(color = "black", linewidth = 1)  # Border for facet labels
  ) +
  facet_wrap(~Region, 
             labeller = labeller(Region = c("North" = "North Taiwan", 
                                            "Green Island" = "Green Island", 
                                            "Xiaoliuqiu" = "Xiaoliuqiu")))  # Modify facet labels



### 7-3. Location & Top 5 Coral OTU---in coral OTU level ####
BC.only <- read.csv("BC.only.csv", header=T, sep=",")

Coral.only.sum <- BC.only %>% 
  group_by(Region, Depth, MajorCategory, OTUs) %>% 
  summarise(Number = length(OTUs)) %>%
  filter(MajorCategory =='Black Corals'| MajorCategory =='Gorgonian Corals'| MajorCategory =='Soft Corals'| MajorCategory =='Stony Corals')

Coral.only.sum$Region <- factor(Coral.only.sum$Region, 
                                levels = c('North', 'Green Island','Xiaoliuqiu'))
Coral.only.sum <- Coral.only.sum[c('Region', 'Depth','OTUs', 'Number')]

Coral.only.percent <- Coral.only.sum %>% group_by(Region, Depth) %>% 
  mutate(Percent = Number/sum(Number))
Coral.only.percent <- Coral.only.percent[c('Region', 'Depth','OTUs', 'Percent')]
# write.csv(Coral.only.percent, 'CoralOTU.percent.csv',row.names = F)

C.top5.otu <- Coral.only.percent %>% slice_max(Percent, n = 5)  # highest 5 coral OTU

# Calculate remaining percentage and add "Other"
C.top5.sum <- C.top5.otu %>% group_by(Region, Depth) %>% summarise(Percent = sum(Percent))
C.other.otu <- C.top5.sum %>% 
  group_by(Region, Depth) %>%
  summarise(Percent = 1 - Percent) %>% 
  mutate(OTUs = "Other")
C.other.otu <- C.other.otu[c('Region', 'Depth','OTUs', 'Percent')]

C.top5.otu <- bind_rows(C.top5.otu, C.other.otu)  # Combine top 5 OTUs with "Other"

C.top5.otu$Region <- factor(C.top5.otu$Region, 
                            levels = c('North', 'Green Island','Xiaoliuqiu'))
C.top5.otu$Depth <- factor(C.top5.otu$Depth, 
                           levels = c('5', '15','30'))
C.top5.otu$OTUs <- factor(C.top5.otu$OTUs, 
                          levels = c('Other','Antipatharia', 'Alcyonacea',
                                     'Xenia/Heteroxenia','Conglomeratusclera','Clavularia','Litophyton','Sinularia/Cladiella/Klyxum',
                                     'Isopora','Anacropora','Montipora','Porites/Montipora','Porites',
                                     'Leptoseris','Pachyseris',
                                     'Turbinaria','Tubastraea','Tubastraea/Cladopsammia','Pocillopora','Echinophyllia/Oxypora',
                                     'Merulinidae.spp1','Merulinidae.spp2','Cyphastrea','Echinopora','Platygyra',
                                     'Millepora'))
# write.csv(C.top5.otu, 'C.top5.coralOTU.csv',row.names = F)


# Top 5 coral plot ####
ggplot(data = C.top5.otu, mapping = aes(x = Depth, y =Percent, fill = OTUs)) +
  geom_bar(position = "stack", stat='identity') + 
  scale_fill_manual(values=c("#edf6f9","#979dac","#BCBFC2",
                             "#FFE3F2","#FFC9E6","#ff9ebb","#ff7aa2","#E5668C",
                             "#ccd5ae","#abdda4","#66c2a5","#9dd9d2","#92c5de",
                             "#d4a373","#DEC3A9",
                             "#72ddf7","#0096c7","#023e8a","#dab6fc","#9381ff",
                             "#fff7bc","#FFDE85","#FFAE52","#FF8952","#FA645F",
                             "#FFE03B")) +
  labs(fill = "Coral OTU") +
  xlab("Depth (m)") +
  scale_y_continuous(name = "Percent (%)",
                     limits = c(0, 1.00),
                     breaks = seq(0, 1.00, 0.2),
                     labels = seq(0, 100, 20)) +
  facet_wrap(~Region, scales = "free_x") +  # Keep x-axis ticks across all facets
  theme(
    text = element_text(size = 12),            # Base text size
    axis.title = element_text(size = 14, face = "bold"),      # Axis titles text size
    axis.text = element_text(size = 12),       # Axis labels text size
    strip.text = element_text(size = 12, face = "bold"),      # Facet labels text size
    legend.text = element_text(size = 6),     # Legend text size
    legend.title = element_text(size = 6),     # Legend title text size
    legend.key.size = unit(0.4, 'cm'),  # Adjusts the overall size of the legend keys
    panel.background = element_blank(),  # Removes the gray background of the plot panel
    plot.background = element_blank(),    # Removes the gray background of the entire plot
    axis.ticks = element_line(color = "black"),  # Keeps the ticks
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),  # Adds black borders around each panel
    strip.background = element_rect(color = "black", linewidth = 1)  # Border for facet labels
  ) +
  facet_wrap(~Region, 
             labeller = labeller(Region = c("North" = "North Taiwan", 
                                            "Green Island" = "Green Island", 
                                            "Xiaoliuqiu" = "Xiaoliuqiu")))  # Modify facet labels



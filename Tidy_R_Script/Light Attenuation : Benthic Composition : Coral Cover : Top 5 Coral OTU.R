### Light Attenuation / Benthic Composition / Top 5 Coral OTU ####
library(ggplot2); library(dplyr); library(scales); library(tidyr)

rm (list = ls())
getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file")

### Light Attenuation ####
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



### Benthic Composition---in Major Category level ####
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



### Location & Top 5 Coral OTU---in coral OTU level ####
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



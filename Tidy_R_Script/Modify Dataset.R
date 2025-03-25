library(dplyr); library(tidyr); library(stringr); library(tidyverse); library(conflicted)
conflict_prefer("select", 'dplyr')
conflict_prefer("filter", "dplyr")

rm (list = ls())
getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/Tidy_file")

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



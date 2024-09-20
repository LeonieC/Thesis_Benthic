library(dplyr); library(tidyr); library(stringr); library(tidyverse); library(conflicted)
conflict_prefer("select", 'dplyr')
conflict_prefer("filter", "dplyr")

getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R")

BC.data <- read.csv("benthic_annotations.csv", header=T, sep=",") %>%
  select(Name, Region, Site, Transect, Depth, Label) %>%
  mutate(Depth = ifelse(Depth == 25, 30, Depth), 
         Location = paste0(Site, '_', Depth))

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
                                'Black coral' = 'Antipatharian',
                                'Soft coral' = 'Octocoral',
                                'Hard coral' = 'Scleractinian',
                                'Actiniaria' = 'Other Sessile Invertebrates',
                                'Ascidian' = 'Other Sessile Invertebrates',
                                'Hydrozoa' = 'Other Sessile Invertebrates',
                                'Zoanthid' = 'Other Sessile Invertebrates',
                                'Corallimorpharia' = 'Other Sessile Invertebrates',
                                'Fish' = 'Other Life',
                                'Other sessile invertebrates' = 'Other Sessile Invertebrates',
                                'Other mobile invertebrates' = 'Other Mobile Invertebrates',
                                'Stable substrate' = 'Stable Substrate',
                                'Unstable substrate' = 'Unstable Substrate', .default = MajorCategory),
         MajorCategory = ifelse(Genus == 'shadow', 'Shadow',
                                ifelse(Genus == 'unknown', 'Unknown', MajorCategory)),
         MajorCategory = recode(MajorCategory, 
                                'Artificial Debris' = 'Other',
                                'Unknown' = 'Other',
                                'TWS' = 'Other', .default = MajorCategory))

BC.fulldata <- BC.data %>% left_join(., BC.labels, by = 'Label')
BC.fulldata <- BC.fulldata[c('Name', 'Region', 'Site', 'Depth', 'Location', 'Transect', 'Label', 'OTUs', 'Genus', 'MajorCategory')]


BC.only <- BC.fulldata %>%   # data only contain the living things
  filter(MajorCategory !='Other' & MajorCategory !='Shadow' & MajorCategory !='Stable Substrate' & MajorCategory !='Unstable Substrate')

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
BC.cover <- BC.cover %>% select(-Image.name, -OTUs, -Morphology, -MorphoFunctionalGroup, -Label)
BC.cover <- BC.cover[c('Region', 'Site', 'Depth', 'Transect', 'Quadrat', 'MajorCategory', 'Genus', 'Cover')]



# output the new create csv
# write.csv(BC.labels, 'BC.labels.csv',row.names = F)
# write.csv(BC.fulldata, 'BC.fulldata.csv',row.names = F)
# write.csv(BC.only, 'BC.only.csv',row.names = F)
# write.csv(BC.cover, 'BC.cover.csv',row.names = F)

# show what Major Category we have in the column: unique(DATA$MajorCategory)
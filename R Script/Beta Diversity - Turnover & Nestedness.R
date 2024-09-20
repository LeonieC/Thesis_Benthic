### Beta Diversity: Turnover & Nestedness of BC & CA ####
## Beta Diversity: Turnover & Nestedness => heatmap---in OTU level
library(dplyr); library(tidyr); library(tibble)
library(betapart); library(graphics)

rm (list = ls())
getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R")

### *all Benthic Community data---in OTU level ####
BC.only <- read.csv("BC.only.csv", header = T, sep = ",")

OTU_bc.sum <- BC.only %>%
  group_by(Region, Location, Depth, MajorCategory, OTUs) %>%
  summarise(Number = length(OTUs))
OTU_bc.sum$Region <- factor(OTU_bc.sum$Region, 
                            levels = c('North', 'Green Island','Xiaoliuqiu'))
OTU_bc.sum$Depth <- factor(OTU_bc.sum$Depth, 
                           levels = c('5', '15','30'))
OTU_bc.sum <- OTU_bc.sum %>%
  mutate(Reg.Dep = factor(paste(Region, Depth, sep = "_"),
                          levels = c("North_5", "North_15", "North_30",
                                     "Green Island_5", "Green Island_15", "Green Island_30",
                                     "Xiaoliuqiu_5", "Xiaoliuqiu_15", "Xiaoliuqiu_30")))

beta.BC <- OTU_bc.sum %>%
  group_by(Reg.Dep, OTUs) %>%
  summarise(Abundance = sum(Number), .groups = 'drop') %>% # remove grouping
  pivot_wider(names_from = OTUs, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  column_to_rownames(var = "Reg.Dep")

# Convert to numeric matrix
beta.BC <- as.matrix(beta.BC)
beta.BC <- beta.BC[levels(OTU_bc.sum$Reg.Dep), ]
# Convert to presence/absence matrix
beta_jar_mat <- ifelse(beta.BC > 1, 1, 0)

# Calculate beta diversity
beta.multi_res <- beta.multi(beta_jar_mat, index.family = "jaccard")
beta_location <- beta.pair(beta_jar_mat, index.family = "jaccard")

# Define a color palette: adjust the colors and the number of levels
my_palette <- colorRampPalette(c("white", "yellow", "orange", "red"))(200)
rev_my_palette <- rev(my_palette)


## Turnover of BC ####
beta_location_tu <- beta_location$beta.jtu
dst <- data.matrix(beta_location_tu)

dst_tu_nor <- dst[1:9, 1:9]  # order => 1:3->N / 4:6->GI / 7:9->XLQ
x_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")


## Nestedness of BC ####
beta_location_ne <- beta_location$beta.jne
dst <- data.matrix(beta_location_ne)

dst_ne_nor <- dst[1:9,1:9] # order => 1:3->N / 4:6->GI / 7:9->XLQ
x_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")


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


## Beta Diversity of BC ####
beta_location_ac <- beta_location$beta.jac
dst <- data.matrix(beta_location_ac)

dst_ac_nor <- dst[1:9,1:9] # order => 1:3->N / 4:6->GI / 7:9->XLQ
x_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")

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



### *only Coral Assemblage data---in OTU level ####
Coral.only <- BC.only %>%
  filter(MajorCategory =='Antipatharian'| MajorCategory =='Gorgonian'| MajorCategory =='Octocoral'| MajorCategory =='Scleractinian')

OTU_coral.sum <- Coral.only %>%
  group_by(Region, Location, Depth, MajorCategory, OTUs) %>%
  summarise(Number = length(OTUs))
OTU_coral.sum$Region <- factor(OTU_coral.sum$Region, 
                               levels = c('North', 'Green Island','Xiaoliuqiu'))
OTU_coral.sum$Depth <- factor(OTU_coral.sum$Depth, 
                              levels = c('5', '15','30'))
OTU_coral.sum <- OTU_coral.sum %>%
  mutate(Reg.Dep = factor(paste(Region, Depth, sep = "_"),
                          levels = c("North_5", "North_15", "North_30",
                                     "Green Island_5", "Green Island_15", "Green Island_30",
                                     "Xiaoliuqiu_5", "Xiaoliuqiu_15", "Xiaoliuqiu_30")))

beta.CA <- OTU_coral.sum %>%
  group_by(Reg.Dep, OTUs) %>%
  summarise(Abundance = sum(Number), .groups = 'drop') %>% # remove grouping
  pivot_wider(names_from = OTUs, values_from = Abundance, values_fill = list(Abundance = 0)) %>%
  column_to_rownames(var = "Reg.Dep")

# Convert to numeric matrix
beta.CA <- as.matrix(beta.CA)
beta.CA <- beta.CA[levels(OTU_coral.sum$Reg.Dep), ]
# Convert to presence/absence matrix
beta_jar_mat.ca <- ifelse(beta.CA > 1, 1, 0)

# Calculate beta diversity
beta.multi_res.ca <- beta.multi(beta_jar_mat.ca, index.family = "jaccard")
beta_location.ca <- beta.pair(beta_jar_mat.ca, index.family = "jaccard")

# Define a color palette: adjust the colors and the number of levels
my_palette <- colorRampPalette(c("white", "yellow", "orange", "red"))(200)
rev_my_palette <- rev(my_palette)


## Turnover of CA ####
beta_location.ca_tu <- beta_location.ca$beta.jtu
dst <- data.matrix(beta_location.ca_tu)

dst_tu_nor <- dst[1:9, 1:9]  # order => 1:3->N / 4:6->GI / 7:9->XLQ
x_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")


## Nestedness of CA ####
beta_location.ca_ne <- beta_location.ca$beta.jne
dst <- data.matrix(beta_location.ca_ne)

dst_ne_nor <- dst[1:9,1:9] # order => 1:3->N / 4:6->GI / 7:9->XLQ
x_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")


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


## Beta Diversity of CA ####
beta_location.ca_ac <- beta_location.ca$beta.jac
dst <- data.matrix(beta_location.ca_ac)

dst_ac_nor <- dst[1:9,1:9] # order => 1:3->N / 4:6->GI / 7:9->XLQ
x_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")
y_labels <- c("N 5", "N 15", "N 30", "GI 5", "GI 15", "GI 30", "XLQ 5", "XLQ 15", "XLQ 30")

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


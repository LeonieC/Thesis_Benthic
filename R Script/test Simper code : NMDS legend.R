### Try code for the indicative species ####

### Simper BC ####
## influential species
simper_res <- simper(OTU_bc.cover, Factor$Region)


# Define the function to extract high contribution species
extract_high_contrib_species <- function(simper_res, threshold = 0.9999) {
  high_contrib_species <- list()
  
  for (comparison in names(simper_res)) {
    species_data <- simper_res[[comparison]]
    
    # Extract cumulative contributions from the species_data
    cumulative_contributions <- species_data$cusum
    
    # Find species above the threshold
    species_above_threshold <- names(cumulative_contributions[cumulative_contributions > threshold])
    high_contrib_species[[comparison]] <- species_above_threshold
  }
  
  return(high_contrib_species)
}

# Apply the function to the simper results
high_contrib_species <- extract_high_contrib_species(simper_res)

# Print the results
print(high_contrib_species)


# NMDS Analysis
nmds1 <- metaMDS(OTU_bc.cover, "bray", type='n')

# Plot NMDS
plot(nmds1, display="sites", type='n', xlab = "NMDS1", ylab = "NMDS2")

# Add points for different regions and depths
points(nmds1, "sites", col = "#e0fbfc", pch = 17, select = Factor$RegionDepth == "North5")
points(nmds1, "sites", col = "#98c1d9", pch = 15, select = Factor$RegionDepth == "North15")
points(nmds1, "sites", col = "#274c77", pch = 16, select = Factor$RegionDepth == "North30")
points(nmds1, "sites", col = "#fff3b0", pch = 17, select = Factor$RegionDepth == "Green Island5")
points(nmds1, "sites", col = "#ffc857", pch = 15, select = Factor$RegionDepth == "Green Island15")
points(nmds1, "sites", col = "#ff9f1c", pch = 16, select = Factor$RegionDepth == "Green Island30")
points(nmds1, "sites", col = "#ecf39e", pch = 17, select = Factor$RegionDepth == "Xiaoliuqiu5")
points(nmds1, "sites", col = "#adc178", pch = 15, select = Factor$RegionDepth == "Xiaoliuqiu15")
points(nmds1, "sites", col = "#709775", pch = 16, select = Factor$RegionDepth == "Xiaoliuqiu30")

# Add hulls for different regions
ordihull(nmds1, groups = Factor$Region, col = c("#ffc857", "#98c1d9", "#adc178"))

# Add stress value
mtext(paste("Stress:", round(nmds1$stress, 4)), side = 1, line = 3, cex = 0.8, adj = 1)

# Add high contribution species
for (comparison in names(high_contrib_species)) {
  species_list <- high_contrib_species[[comparison]]
  
  # Extract species scores for NMDS plot
  species_scores <- scores(nmds1, display = "species")
  
  # Filter species that are in the list of high contribution species
  filtered_scores <- species_scores[rownames(species_scores) %in% species_list, ]
  
  # Add species to plot
  text(filtered_scores[, 1], filtered_scores[, 2], labels = rownames(filtered_scores), cex = 0.5, pos = 4, col = "#979dac")
}



### Simper CA ####
## influential species
simper_res.ca <- simper(OTU_coral.cover, Factor2$Region)


# Define the function to extract high contribution species
extract_high_contrib_species.ca <- function(simper_res.ca, threshold = 0.9999) {
  high_contrib_species.ca <- list()
  
  for (comparison in names(simper_res.ca)) {
    species_data.ca <- simper_res.ca[[comparison]]
    
    # Extract cumulative contributions from the species_data.ca
    cumulative_contributions.ca <- species_data.ca$cusum
    
    # Find species above the threshold
    species_above_threshold.ca <- names(cumulative_contributions.ca[cumulative_contributions.ca > threshold])
    high_contrib_species.ca[[comparison]] <- species_above_threshold.ca
  }
  
  return(high_contrib_species.ca)
}

# Apply the function to the simper results
high_contrib_species.ca <- extract_high_contrib_species.ca(simper_res.ca)

# Print the results
print(high_contrib_species.ca)


# NMDS Analysis
nmds2 <- metaMDS(OTU_coral.cover, "bray", type='n')

# Plot NMDS
plot(nmds2, display="sites", type='n', xlab = "NMDS1", ylab = "NMDS2")

# Add points for different regions and depths
points(nmds2, "sites", col = "#e0fbfc", pch = 17, select = Factor$RegionDepth == "North5")
points(nmds2, "sites", col = "#98c1d9", pch = 15, select = Factor$RegionDepth == "North15")
points(nmds2, "sites", col = "#274c77", pch = 16, select = Factor$RegionDepth == "North30")
points(nmds2, "sites", col = "#fff3b0", pch = 17, select = Factor$RegionDepth == "Green Island5")
points(nmds2, "sites", col = "#ffc857", pch = 15, select = Factor$RegionDepth == "Green Island15")
points(nmds2, "sites", col = "#ff9f1c", pch = 16, select = Factor$RegionDepth == "Green Island30")
points(nmds2, "sites", col = "#ecf39e", pch = 17, select = Factor$RegionDepth == "Xiaoliuqiu5")
points(nmds2, "sites", col = "#adc178", pch = 15, select = Factor$RegionDepth == "Xiaoliuqiu15")
points(nmds2, "sites", col = "#709775", pch = 16, select = Factor$RegionDepth == "Xiaoliuqiu30")

# Add hulls for different regions
ordihull(nmds2, groups = Factor$Region, col = c("#ffc857", "#98c1d9", "#adc178"))

# Add stress value
mtext(paste("Stress:", round(nmds2$stress, 4)), side = 1, line = 3, cex = 0.8, adj = 1)

# Add high contribution species
for (comparison in names(high_contrib_species.ca)) {
  species_list <- high_contrib_species.ca[[comparison]]
  
  # Extract species scores for NMDS plot
  species_scores <- scores(nmds2, display = "species")
  
  # Filter species that are in the list of high contribution species
  filtered_scores <- species_scores[rownames(species_scores) %in% species_list, ]
  
  # Add species to plot
  text(filtered_scores[, 1], filtered_scores[, 2], labels = rownames(filtered_scores), cex = 0.5, pos = 4, col = "#979dac")
}



### Legend for NMDS plot ####
library(fields)
par(mar = c(8, 0, 2, 0)) # Adjust margins
image.plot()
# Legend for Depth
legend("topright", 
       legend = c("5", "15", "30"),
       pch = c(17, 15, 16),
       title = "Depth")

legend("bottomleft", 
       legend = c("North", "Green Island", "Xiaoliuqiu"),
       col = c("#98c1d9", "#ffc857", "#adc178"),
       title = "Region",
       ncol = 3)



# Add a legend for Depth (symbols)
legend("topright", legend = c("5", "15", "30"), pch = c(17, 15, 16), col = "black", title = "Depth")


# Add a legend for Region (colors)
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 1)))

# Create legend items
legend_items <- c("North", "Green Island", "Xiaoliuqiu")
colors <- c("#98c1d9", "#ffc857", "#adc178")
text_colors <- c("#98c1d9", "#ffc857", "#adc178")

# Draw legend
grid.rect(x = 0.5, y = 0.5, width = 0.9, height = 0.9, gp = gpar(fill = "white"))
for (i in 1:length(legend_items)) {
  grid.text(legend_items[i], x = 0.1, y = 0.9 - i*0.1, gp = gpar(col = text_colors[i]), just = "left")
  grid.points(x = 0.02, y = 0.9 - i*0.1, gp = gpar(col = colors[i]))
}

# Draw legend background
grid.rect(x = 0.5, y = 0.5, width = 0.9, height = 0.9, gp = gpar(fill = "white"))

# Draw the title
grid.text("Region", x = 0.2, y = 0.85, gp = gpar(fontsize = 14, fontface = "bold"))

# Draw legend items
for (i in 1:length(legend_items)) {
  grid.text(legend_items[i], x = 0.1, y = 0.8 - i*0.1, gp = gpar(col = text_colors[i], fontface = "bold"), just = "left")
  grid.points(x = 0.02, y = 0.9 - i*0.1, gp = gpar(col = colors[i]))
}

### Map for Study Sites ####
library(sf); library(ggplot2); library(maps); library(ggspatial)
library(dplyr); library(tidyr); library(stringr); library(tidyverse)

getwd()
setwd("/Users/leonie/Desktop/R_for_Benthic/Benthic_R")

# North
#FL 25.04, 121.94
#BT 25.13, 121.85
#KI 25.19, 121.78

# Green Island
#GW 22.64, 121.48
#DBS 22.63, 121.49
#GG 22.67, 121.49

# Xiaoliuqiu
#LC 22.34, 120.38
#VN 22.33, 120.35
#DF 22.32, 120.37



## Taiwan ####
# Load the shapefile
taiwan_shapefile <- st_read("/Users/leonie/Desktop/R_for_Benthic/Benthic_R/gadm41_TWN_0.shp")  # Adjust the path and file name accordingly

# Create a data frame with longitude and latitude coordinates
points_df <- data.frame(
  Longitude = c(121.95, 121.85 , 121.78,
                121.48, 121.49, 121.49,
                120.392, 120.355, 120.375),  # longitudes
  Latitude = c(25.04, 25.14, 25.19,
               22.64, 22.635, 22.684,
               22.345, 22.34, 22.33),  # latitudes
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
  geom_rect(aes(xmin = 120.28, xmax = 120.45, ymin = 22.25, ymax = 22.40), color = "#adc178", fill = NA, linetype = "solid") +  # Xiaoliuqiu
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

# Convert to sf object
north_sf <- st_as_sf(north_df, coords = c("Longitude", "Latitude"), crs = 4326)

# Plot the inset map for North
inset_north <- ggplot(data = north_sf) +
  geom_sf(data = taiwan_shapefile, fill = "lightgray", color = "black") +
  geom_sf(fill = "blue", color = "blue", size = 2) +
  coord_sf(xlim = c(121.68, 122.08), ylim = c(24.95, 25.25)) +
  geom_text(data = north_df, aes(x = Longitude, y = Latitude, label = site), vjust = -1, color = "blue") +
  theme_void()

print(inset_north)



## Green Island ####
# Create a data frame for Green Island
green_island_df <- points_df %>% filter(region == "Green Island")

# Convert to sf object
green_island_sf <- st_as_sf(green_island_df, coords = c("Longitude", "Latitude"), crs = 4326)

# Plot the inset map for Green Island
inset_green_island <- ggplot(data = green_island_sf) +
  geom_sf(data = taiwan_shapefile, fill = "lightgray", color = "black") +
  geom_sf(fill = "blue", color = "blue", size = 2) +
  coord_sf(xlim = c(121.45, 121.53), ylim = c(22.60, 22.70)) +
  geom_text(data = green_island_df, aes(x = Longitude, y = Latitude, label = site), vjust = -1, color = "blue") +
  theme_void()

print(inset_green_island)



## Xiaoliuqiu ####
# Create a data frame for Xiaoliuqiu
xiaoliuqiu_df <- points_df %>% filter(region == "Xiaoliuqiu")

# Convert to sf object
xiaoliuqiu_sf <- st_as_sf(xiaoliuqiu_df, coords = c("Longitude", "Latitude"), crs = 4326)

# Plot the inset map for Xiaoliuqiu
inset_xiaoliuqiu <- ggplot(data = xiaoliuqiu_sf) +
  geom_sf(data = taiwan_shapefile, fill = "lightgray", color = "black") +
  geom_sf(fill = "blue", color = "blue", size = 2) +
  coord_sf(xlim = c(120.33, 120.40), ylim = c(22.30, 22.38)) +
  geom_text(data = xiaoliuqiu_df, aes(x = Longitude, y = Latitude, label = site), vjust = -1, color = "blue") +
  theme_void()

print(inset_xiaoliuqiu)



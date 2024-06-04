# New script for PAME regions image----

# Loading essential packages----
library(sf)
library(sp)

# Loading data----

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data"
mylocs <- readRDS(paste0(datadir, "/SEATRACK_FUGLA_20220307_v2.3_FA.rds"))
summary_info <- readRDS(paste0(datadir,"/summaryTable.rds"))
indiv_merged_df <- merge(mylocs, summary_info, by = "ring") 
nbs_mylocs <- indiv_merged_df %>% dplyr::mutate(year = year(timestamp))
names(nbs_mylocs)[names(nbs_mylocs) == "ring"] <- "individ_id"
all_data <- nbs_mylocs
df <- all_data[!is.na(all_data$timestamp), ]
df <- df %>% dplyr::mutate(month = month(timestamp))
months <- sort(unique(df$month))
df_mod <- df %>% mutate(tracking_year = ifelse(month < 7, year - 1, year))
df_mod_sf <- df_mod %>% 
  sf::st_as_sf(coords = c("lon", "lat"), crs = 4326) 
ind_merge_sf <- df_mod_sf

datadir <- "/Users/ameydanole/Desktop/ENS_Rennes/argh/Microplastic_ingestion_by_fulmarus_glacialis/Ind/input_data"
PAME_shapefile <- st_read(paste0(datadir,"/PAME/modified_LME.shp"))
PAME_curated <- PAME_shapefile[c(2,8,9,12,14),]
PAME_curated_without_Faroe <- PAME_curated[-2,]

worldmap <- ggplot2::map_data('world')
water_polygons <- rnaturalearth::ne_download(
  category = "physical",
  type = "ocean",
  scale = 110,
  returnclass = "sf"
)

# Plots----

PAME_plot <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill = "lightgrey") +
  geom_sf(data = water_polygons, fill = "lightblue", color = "black") +  # Overlay water polygons
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Set bold face for plot title
    axis.title = element_text(face = "bold"),  # Set bold face for axis titles
    legend.title = element_text(face = "bold"),  # Set bold face for legend title
    legend.text = element_text(face = "bold")  # Set bold face for legend text
  ) +
  ggtitle("Relevant Protection of the Arctic Marine Environment (PAME) zones") +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_sf(data = PAME_shapefile, aes(), fill = NA, color = "black", size = 7) +
  geom_sf(data = PAME_curated, aes(), fill = "white", color = "black", size = 7) +
  geom_sf_label(data = PAME_curated, aes(label = LME_NAME), 
                size = 3.5, color = "black", fill = "white", label.padding = unit(0.1, "lines"),
                label.size = 0.5) +  # Adjust size and aesthetics of the label
  geom_point(data = ind_merge_sf, aes(x = col_lon, y = col_lat, color = colony)) + 
  coord_sf(xlim = c(-95, 95), ylim = c(30, 90)) +  # Set limits of the plot
  labs(color = "Colony")  # Set legend title
print(PAME_plot)

# Dealing with Faroe Plateau---- 

PAME_plot <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill = "lightgrey") +
  geom_sf(data = water_polygons, fill = "lightblue", color = "black") +  # Overlay water polygons
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Set bold face for plot title
    axis.title = element_text(face = "bold"),  # Set bold face for axis titles
    legend.title = element_text(face = "bold"),  # Set bold face for legend title
    legend.text = element_text(face = "bold")  # Set bold face for legend text
  ) +
  ggtitle("Relevant Protection of the Arctic Marine Environment (PAME) zones") +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_sf(data = PAME_shapefile, aes(), fill = NA, color = "black", size = 7) +
  geom_sf(data = PAME_curated, aes(), fill = "white", color = "black", size = 7) +
  geom_label(data = PAME_curated %>% filter(LME_NAME == "Faroe Plateau"),
             aes(x = -8, y = 63, label = LME_NAME), 
             size = 3.5, color = "black",
             label.padding = unit(0.1, "lines"), label.size = 0.5,
             nudge_x = -10, nudge_y = -2) +  # Adjust label position for Faroe Islands
  geom_sf_label(data = PAME_curated_without_Faroe, aes(label = LME_NAME), 
                size = 3.5, color = "black", fill = "white", label.padding = unit(0.1, "lines"),
                label.size = 0.5) +  # Adjust size and aesthetics of the label
  geom_point(data = ind_merge_sf, aes(x = col_lon, y = col_lat, color = colony)) + 
  coord_sf(xlim = c(-95, 95), ylim = c(30, 90)) +  # Set limits of the plot
  labs(color = "Colony")  # Set legend title
print(PAME_plot)

# Now final plot with colony colours accorded----

names(ind_merge_sf)[5] <- "Colony"
to_combine <- c("Skjalfandi", "Langanes")
ind_merge_sf$Colony[ind_merge_sf$Colony %in% to_combine] <- "Iceland"

# Create a data frame with the colony colors
colony_colors <- data.frame(
  Colony = c("Little Saltee", "Inishkea", "Isle of Canna", "Eynhallow", "Jarsteinen",
             "Faroe Islands", "Iceland","Jan Mayen", "Bjørnøya", "Alkefjellet"),
  Color = c("#440154FF", "#482878FF", "#3E4A89FF", "#31688EFF", "#26828EFF",
            "#1F9E89FF", "#35B779FF", "#6DCD59FF", "#B4DE2CFF", "#FDE725FF")
)


# Merge colony colors with ind_merge_sf data
ind_merge_sf <- merge(ind_merge_sf, colony_colors, by = "Colony", all.x = TRUE)
ind_merge_sf$Colony <- as.factor(ind_merge_sf$Colony)
order_to_match <- c("Alkefjellet", "Bjørnøya", "Jan Mayen", "Iceland", "Faroe Islands", "Jarsteinen", "Eynhallow", "Isle of Canna", "Inishkea", "Little Saltee")
ind_merge_sf$Colony <- factor(ind_merge_sf$Colony, levels = order_to_match)

# Final plot
  PAME_plot <- ggplot() +
  geom_polygon(data = worldmap, aes(x = long, y = lat, group = group), fill = "lightgrey") +
  geom_sf(data = water_polygons, fill = "lightblue", color = "black") +  # Overlay water polygons
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"), # Set bold face for plot title
    axis.title = element_text(face = "bold"),  # Set bold face for axis titles
    legend.title = element_text(face = "bold"),  # Set bold face for legend title
    legend.text = element_text(face = "bold")  # Set bold face for legend text
  ) +
  ggtitle("Protection of the Arctic Marine Environment (PAME) zones and corresponding EcoQO values") +
  xlab("Longitude") +
  ylab("Latitude") +
  geom_sf(data = PAME_shapefile, aes(), fill = NA, color = "black", size = 7) +
  geom_sf(data = PAME_curated, aes(), fill = "white", color = "black", size = 7) +
  geom_label(data = PAME_curated %>% filter(LME_NAME == "Faroe Plateau"),
             aes(x = -8, y = 63, label = LME_NAME), 
             size = 3.5, color = "black",
             label.padding = unit(0.1, "lines"), label.size = 0.5,
             nudge_x = -10, nudge_y = -2) +  # Adjust label position for Faroe Islands
  geom_sf_label(data = PAME_curated_without_Faroe, aes(label = LME_NAME), 
                size = 3.5, color = "black", fill = "white", label.padding = unit(0.1, "lines"),
                label.size = 0.5) +  # Adjust size and aesthetics of the label
  geom_point(data = ind_merge_sf, aes(x = col_lon, y = col_lat), shape = 17, size = 6, color = "black") +  # Overlay larger black triangles with a smaller size
  geom_point(data = ind_merge_sf, aes(x = col_lon, y = col_lat, color = Colony), shape = 17, size = 5) + 
  scale_color_manual(values = rev(colony_colors$Color)) + 
  coord_sf(xlim = c(-95, 95), ylim = c(30, 90)) +  # Set limits of the plot
  labs(color = "Colony")  # Set legend title
print(PAME_plot)
ggsave("Relevant_PAME_zones_with_colonies_and_EcoQO.png", height = 8, width= 12, unit = "in", dpi = 900, plot = PAME_plot, bg = "white")



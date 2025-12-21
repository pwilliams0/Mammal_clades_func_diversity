#===================================#
# Maps of global diversity patterns #
#===================================#

library(tidyverse)
library(sf)
library(terra)
library(colorspace)
library(rnaturalearth)
library(rnaturalearthdata)

# ----- LOAD DATA ----------

# Load grid cells
cells <- st_read("Data/cells.shp")
# Load world outline
world <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_transform(crs=crs(cells))

# Load diversity data
rodent_div_metrics <- read.csv("Data/rodent_div_metrics.csv")
chirop_div_metrics <- read.csv("Data/chirop_div_metrics.csv")
eulipo_div_metrics <- read.csv("Data/eulipo_div_metrics.csv")
primat_div_metrics <- read.csv("Data/primat_div_metrics.csv")
carniv_div_metrics <- read.csv("Data/carniv_div_metrics.csv")
artiod_div_metrics <- read.csv("Data/artiod_div_metrics.csv")

# ----- SPECIES RICHNESS & FUNCTIONAL DIVERSITY ----------

# Set color palette
colorspace::sequential_hcl(n = 7, h = 12, c = c(50, 80, NA),
                           l = c(20, 97), power = 1,
                           register = "palette")

# Plot alpha diversity maps
# Arguments:
#   div_data = table with diversity of each grid cell
#   column = diversity metric column name
#   cells = table of grid cells with cell_id and geometry
#   palette = color palette to use
diversity_map <- function(div_data, column, cells) {
  # Format data
  data <- div_data %>%
    left_join(cells, "cell_id") %>%
    mutate(div = get(column)) %>%
    dplyr::filter(div > -1000)
  # Plot map of diversity metric
  map <- ggplot() +
    geom_sf(data=world, colour="grey85", linewidth=.1) +
    geom_sf(data=data, linewidth=0.05,
            aes(fill=div, color=div, geometry=geometry)) +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.position = "bottom") +
    theme(legend.key.height = unit(.12, "cm"),
          legend.key.width = unit(.8,"cm"))
    
  if(min(data$div) < 0){
    map <- map + scale_fill_continuous_diverging("Blue-Red 3") +
      scale_color_continuous_diverging("Blue-Red 3")} else {
        map <- map + scale_fill_continuous_sequential("Palette") +
          scale_color_continuous_sequential("Palette")}
  map
}

# Rodentia ----------
diversity_map(rodent_div_metrics, "SRic", cells)
ggsave("Results/DiversityMaps/map_rodent_SRic.svg",
       width = 2, height = 1.2)
diversity_map(rodent_div_metrics, "SES.FRic", cells)
ggsave("Results/DiversityMaps/map_rodent_SES.FRic.svg",
       width = 2, height = 1.2)
diversity_map(rodent_div_metrics, "FDis", cells)
ggsave("Results/DiversityMaps/map_rodent_FDis.svg",
       width = 2, height = 1.2)

# Chiroptera ----------
diversity_map(chirop_div_metrics, "SRic", cells)
ggsave("Results/DiversityMaps/map_chirop_SRic.svg",
       width = 2, height = 1.2)
diversity_map(chirop_div_metrics, "SES.FRic", cells)
ggsave("Results/DiversityMaps/map_chirop_SES.FRic.svg",
       width = 2, height = 1.2)
diversity_map(chirop_div_metrics, "FDis", cells)
ggsave("Results/DiversityMaps/map_chirop_FDis.svg",
       width = 2, height = 1.2)

# Eulipotyphla ----------
diversity_map(eulipo_div_metrics, "SRic", cells)
ggsave("Results/DiversityMaps/map_eulipo_SRic.svg",
       width = 2, height = 1.2)
diversity_map(eulipo_div_metrics, "SES.FRic", cells)
ggsave("Results/DiversityMaps/map_eulipo_SES.FRic.svg",
       width = 2, height = 1.2)
diversity_map(eulipo_div_metrics, "FDis", cells)
ggsave("Results/DiversityMaps/map_eulipo_FDis.svg",
       width = 2, height = 1.2)

# Primates ----------
diversity_map(primat_div_metrics, "SRic", cells)
ggsave("Results/DiversityMaps/map_primat_SRic.svg",
       width = 2, height = 1.2)
# Manually set limits for SES FRic
data <- primat_div_metrics %>%
  left_join(cells, "cell_id") %>%
  mutate(div = SES.FRic) %>%
  dplyr::filter(div > -1000)
map <- ggplot() +
  geom_sf(data=world, colour="grey85", linewidth=.1) +
  geom_sf(data=data, linewidth=0.05,
          aes(fill=div, color=div, geometry=geometry)) +
  theme_void() +
  theme(legend.title = element_blank(),
        legend.text = element_text(size=7),
        legend.position = "bottom") +
  theme(legend.key.height = unit(.12, "cm"),
        legend.key.width = unit(.7,"cm")) +
  scale_fill_continuous_diverging("Blue-Red 3",
                                  limits = c(min(data$div), 7)) +
  scale_color_continuous_diverging("Blue-Red 3",
                                   limits = c(min(data$div), 7))
map
ggsave("Results/DiversityMaps/map_primat_SES.FRic.svg",
       width = 2, height = 1.2)
diversity_map(primat_div_metrics, "FDis", cells)
ggsave("Results/DiversityMaps/map_primat_FDis.svg",
       width = 2, height = 1.2)

# Artiodactyla ----------
diversity_map(artiod_div_metrics, "SRic", cells)
ggsave("Results/DiversityMaps/map_artiod_SRic.svg",
       width = 2, height = 1.2)
diversity_map(artiod_div_metrics, "SES.FRic", cells)
ggsave("Results/DiversityMaps/map_artiod_SES.FRic.svg",
       width = 2, height = 1.2)
diversity_map(artiod_div_metrics, "FDis", cells)
ggsave("Results/DiversityMaps/map_artiod_FDis.svg",
       width = 2, height = 1.2)

# Carnivora ----------
diversity_map(carniv_div_metrics, "SRic", cells)
ggsave("Results/DiversityMaps/map_carniv_SRic.svg",
       width = 2, height = 1.2)
diversity_map(carniv_div_metrics, "SES.FRic", cells)
ggsave("Results/DiversityMaps/map_carniv_SES.FRic.svg",
       width = 2, height = 1.2)
diversity_map(carniv_div_metrics, "FDis", cells)
ggsave("Results/DiversityMaps/map_carniv_FDis.svg",
       width = 2, height = 1.2)

# ----- TRAIT DISTRIBUTIONS ----------

# Calculate mean and sd of a given trait for each cell
# Arguments:
#   trait = trait to calculate mean and sd
#   trait_df = dataframe containing species trait data
#   comm_mat = community matrix with presence/absence of species in each cell
#   div_metrics = diversity metrics of each cell
calc_trait <- function(trait, trait_df, comm_mat, div_metrics) {
  data <- comm_mat %>%
    rownames_to_column("cell_id") %>%
    # Create separate rows for each combination of species and cell
    gather(key = "Binomial.1.2", value = "present", -cell_id) %>%
    filter(present == 1) %>%
    # Remove cells that had too few species to calculate functional diversity
    right_join(div_metrics, "cell_id") %>%
    filter(FDis >= 0) %>%
    # Join with species traits
    left_join(trait_df, by = "Binomial.1.2") %>%
    mutate(logMass = log(Mass.g)) %>%
    # Calculate mean and sd of the trait for each cell
    group_by(cell_id) %>%
    mutate(mean_trait = mean(get(trait)),
           sd_trait = sd(get(trait))) %>%
    slice(1) %>%
    select(cell_id, mean_trait, sd_trait) %>%
    # Join with "cells" to get geometry
    left_join(cells, "cell_id")
  
  data
}

# Map mean or sd trait values of each cell
# Arguments:
#   data = dataframe containing mean and sd trait values for each cell
#   val = value to plot, either "mean_trait" or "sd_trait"
map_trait <- function(data, val) {
  ggplot() +
    geom_sf(data=world, colour="grey85", linewidth=.1) +
    geom_sf(data=data, linewidth=0.05,
            aes(fill=get(val), color=get(val), geometry=geometry)) +
    theme_void() +
    scale_fill_viridis_c(option = "magma", direction = -1) +
    scale_color_viridis_c(option = "magma", direction = -1) +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.position = "bottom") +
    theme(legend.key.height = unit(.12, "cm"),
          legend.key.width = unit(.8,"cm"))
}

# Rodentia ----------
trait_df <- read.csv("Data/traits_rodent.csv")
comm_mat <- readRDS("Data/comm_mat_rodent.RDS")
div_metrics <- rodent_div_metrics
# logMass
logMass <- calc_trait("logMass", trait_df, comm_mat, div_metrics)
map_trait(logMass, "mean_trait")
ggsave("Results/TraitMaps/rodent_logMass_mean.svg",
       width = 2, height = 1.2)
map_trait(logMass, "sd_trait")
ggsave("Results/TraitMaps/rodent_logMass_sd.svg",
       width = 2, height = 1.2)
# ForStrat.Arboreal
Arboreal <- calc_trait("ForStrat.Arboreal", trait_df, comm_mat, div_metrics)
map_trait(Arboreal, "mean_trait")
ggsave("Results/TraitMaps/rodent_Arboreal_mean.svg",
       width = 2, height = 1.2)
map_trait(Arboreal, "sd_trait")
ggsave("Results/TraitMaps/rodent_Arboreal_sd.svg",
       width = 2, height = 1.2)
# ForStrat.Ground
Ground <- calc_trait("ForStrat.Ground", trait_df, comm_mat, div_metrics)
map_trait(Ground, "mean_trait")
ggsave("Results/TraitMaps/rodent_Ground_mean.svg",
       width = 2, height = 1.2)
map_trait(Ground, "sd_trait")
ggsave("Results/TraitMaps/rodent_Ground_sd.svg",
       width = 2, height = 1.2)
# Diet.Invertebrate
Invert <- calc_trait("Diet.Invertebrate", trait_df, comm_mat, div_metrics)
map_trait(Invert, "mean_trait")
ggsave("Results/TraitMaps/rodent_Invert_mean.svg",
       width = 2, height = 1.2)
map_trait(Invert, "sd_trait")
ggsave("Results/TraitMaps/rodent_Invert_sd.svg",
       width = 2, height = 1.2)
# Diet.Vertebrate
Vert <- calc_trait("Diet.Vertebrate", trait_df, comm_mat, div_metrics)
map_trait(Vert, "mean_trait")
ggsave("Results/TraitMaps/rodent_Vert_mean.svg",
       width = 2, height = 1.2)
map_trait(Vert, "sd_trait")
ggsave("Results/TraitMaps/rodent_Vert_sd.svg",
       width = 2, height = 1.2)
# Diet.Fruit
Fruit <- calc_trait("Diet.Fruit", trait_df, comm_mat, div_metrics)
map_trait(Fruit, "mean_trait")
ggsave("Results/TraitMaps/rodent_Fruit_mean.svg",
       width = 2, height = 1.2)
map_trait(Fruit, "sd_trait")
ggsave("Results/TraitMaps/rodent_Fruit_sd.svg",
       width = 2, height = 1.2)
# Diet.Seed
Seed <- calc_trait("Diet.Seed", trait_df, comm_mat, div_metrics)
map_trait(Seed, "mean_trait")
ggsave("Results/TraitMaps/rodent_Seed_mean.svg",
       width = 2, height = 1.2)
map_trait(Seed, "sd_trait")
ggsave("Results/TraitMaps/rodent_Seed_sd.svg",
       width = 2, height = 1.2)
# Diet.Herb
Herb <- calc_trait("Diet.Herb", trait_df, comm_mat, div_metrics)
map_trait(Herb, "mean_trait")
ggsave("Results/TraitMaps/rodent_Herb_mean.svg",
       width = 2, height = 1.2)
map_trait(Herb, "sd_trait")
ggsave("Results/TraitMaps/rodent_Herb_sd.svg",
       width = 2, height = 1.2)

# Chiroptera ----------
trait_df <- read.csv("Data/traits_chirop.csv")
comm_mat <- readRDS("Data/comm_mat_chirop.RDS")
div_metrics <- chirop_div_metrics
# logMass
logMass <- calc_trait("logMass", trait_df, comm_mat, div_metrics)
map_trait(logMass, "mean_trait")
ggsave("Results/TraitMaps/chirop_logMass_mean.svg",
       width = 2, height = 1.2)
map_trait(logMass, "sd_trait")
ggsave("Results/TraitMaps/chirop_logMass_sd.svg",
       width = 2, height = 1.2)
# ForStrat.Arboreal
Arboreal <- calc_trait("ForStrat.Arboreal", trait_df, comm_mat, div_metrics)
map_trait(Arboreal, "mean_trait")
ggsave("Results/TraitMaps/chirop_Arboreal_mean.svg",
       width = 2, height = 1.2)
map_trait(Arboreal, "sd_trait")
ggsave("Results/TraitMaps/chirop_Arboreal_sd.svg",
       width = 2, height = 1.2)
# ForStrat.Ground
Ground <- calc_trait("ForStrat.Ground", trait_df, comm_mat, div_metrics)
map_trait(Ground, "mean_trait")
ggsave("Results/TraitMaps/chirop_Ground_mean.svg",
       width = 2, height = 1.2)
map_trait(Ground, "sd_trait")
ggsave("Results/TraitMaps/chirop_Ground_sd.svg",
       width = 2, height = 1.2)
# ForStrat.Aerial
Aerial <- calc_trait("ForStrat.Aerial", trait_df, comm_mat, div_metrics)
map_trait(Aerial, "mean_trait")
ggsave("Results/TraitMaps/chirop_Aerial_mean.svg",
       width = 2, height = 1.2)
map_trait(Aerial, "sd_trait")
ggsave("Results/TraitMaps/chirop_Aerial_sd.svg",
       width = 2, height = 1.2)
# Diet.Invertebrate
Invert <- calc_trait("Diet.Invertebrate", trait_df, comm_mat, div_metrics)
map_trait(Invert, "mean_trait")
ggsave("Results/TraitMaps/chirop_Invert_mean.svg",
       width = 2, height = 1.2)
map_trait(Invert, "sd_trait")
ggsave("Results/TraitMaps/chirop_Invert_sd.svg",
       width = 2, height = 1.2)
# Diet.Vertebrate
Vert <- calc_trait("Diet.Vertebrate", trait_df, comm_mat, div_metrics)
map_trait(Vert, "mean_trait")
ggsave("Results/TraitMaps/chirop_Vert_mean.svg",
       width = 2, height = 1.2)
map_trait(Vert, "sd_trait")
ggsave("Results/TraitMaps/chirop_Vert_sd.svg",
       width = 2, height = 1.2)
# Diet.Fruit
Fruit <- calc_trait("Diet.Fruit", trait_df, comm_mat, div_metrics)
map_trait(Fruit, "mean_trait")
ggsave("Results/TraitMaps/chirop_Fruit_mean.svg",
       width = 2, height = 1.2)
map_trait(Fruit, "sd_trait")
ggsave("Results/TraitMaps/chirop_Fruit_sd.svg",
       width = 2, height = 1.2)
# Diet.Nect
Nect <- calc_trait("Diet.Nect", trait_df, comm_mat, div_metrics)
map_trait(Nect, "mean_trait")
ggsave("Results/TraitMaps/chirop_Nect_mean.svg",
       width = 2, height = 1.2)
map_trait(Nect, "sd_trait")
ggsave("Results/TraitMaps/chirop_Nect_sd.svg",
       width = 2, height = 1.2)

# Eulipotyphla ----------
trait_df <- read.csv("Data/traits_eulipo.csv")
comm_mat <- readRDS("Data/comm_mat_eulipo.RDS")
div_metrics <- eulipo_div_metrics
# logMass
logMass <- calc_trait("logMass", trait_df, comm_mat, div_metrics)
map_trait(logMass, "mean_trait")
ggsave("Results/TraitMaps/eulipo_logMass_mean.svg",
       width = 2, height = 1.2)
map_trait(logMass, "sd_trait")
ggsave("Results/TraitMaps/eulipo_logMass_sd.svg",
       width = 2, height = 1.2)
# ForStrat.Arboreal
Arboreal <- calc_trait("ForStrat.Arboreal", trait_df, comm_mat, div_metrics)
map_trait(Arboreal, "mean_trait")
ggsave("Results/TraitMaps/eulipo_Arboreal_mean.svg",
       width = 2, height = 1.2)
map_trait(Arboreal, "sd_trait")
ggsave("Results/TraitMaps/eulipo_Arboreal_sd.svg",
       width = 2, height = 1.2)
# ForStrat.Ground
Ground <- calc_trait("ForStrat.Ground", trait_df, comm_mat, div_metrics)
map_trait(Ground, "mean_trait")
ggsave("Results/TraitMaps/eulipo_Ground_mean.svg",
       width = 2, height = 1.2)
map_trait(Ground, "sd_trait")
ggsave("Results/TraitMaps/eulipo_Ground_sd.svg",
       width = 2, height = 1.2)
# Diet.Invertebrate
Invert <- calc_trait("Diet.Invertebrate", trait_df, comm_mat, div_metrics)
map_trait(Invert, "mean_trait")
ggsave("Results/TraitMaps/eulipo_Invert_mean.svg",
       width = 2, height = 1.2)
map_trait(Invert, "sd_trait")
ggsave("Results/TraitMaps/eulipo_Invert_sd.svg",
       width = 2, height = 1.2)
# Diet.Vertebrate
Vert <- calc_trait("Diet.Vertebrate", trait_df, comm_mat, div_metrics)
map_trait(Vert, "mean_trait")
ggsave("Results/TraitMaps/eulipo_Vert_mean.svg",
       width = 2, height = 1.2)
map_trait(Vert, "sd_trait")
ggsave("Results/TraitMaps/eulipo_Vert_sd.svg",
       width = 2, height = 1.2)
# Diet.Plant
Plant <- calc_trait("Diet.Plant", trait_df, comm_mat, div_metrics)
map_trait(Plant, "mean_trait")
ggsave("Results/TraitMaps/eulipo_Plant_mean.svg",
       width = 2, height = 1.2)
map_trait(Plant, "sd_trait")
ggsave("Results/TraitMaps/eulipo_Plant_sd.svg",
       width = 2, height = 1.2)

# Primates ----------
trait_df <- read.csv("Data/traits_primat.csv")
comm_mat <- readRDS("Data/comm_mat_primat.RDS")
div_metrics <- primat_div_metrics
# logMass
logMass <- calc_trait("logMass", trait_df, comm_mat, div_metrics)
map_trait(logMass, "mean_trait")
ggsave("Results/TraitMaps/primat_logMass_mean.svg",
       width = 2, height = 1.2)
map_trait(logMass, "sd_trait")
ggsave("Results/TraitMaps/primat_logMass_sd.svg",
       width = 2, height = 1.2)
# ForStrat.Arboreal
Arboreal <- calc_trait("ForStrat.Arboreal", trait_df, comm_mat, div_metrics)
map_trait(Arboreal, "mean_trait")
ggsave("Results/TraitMaps/primat_Arboreal_mean.svg",
       width = 2, height = 1.2)
map_trait(Arboreal, "sd_trait")
ggsave("Results/TraitMaps/primat_Arboreal_sd.svg",
       width = 2, height = 1.2)
# ForStrat.Ground
Ground <- calc_trait("ForStrat.Ground", trait_df, comm_mat, div_metrics)
map_trait(Ground, "mean_trait")
ggsave("Results/TraitMaps/primat_Ground_mean.svg",
       width = 2, height = 1.2)
map_trait(Ground, "sd_trait")
ggsave("Results/TraitMaps/primat_Ground_sd.svg",
       width = 2, height = 1.2)
# Diet.Invertebrate
Invert <- calc_trait("Diet.Invertebrate", trait_df, comm_mat, div_metrics)
map_trait(Invert, "mean_trait")
ggsave("Results/TraitMaps/primat_Invert_mean.svg",
       width = 2, height = 1.2)
map_trait(Invert, "sd_trait")
ggsave("Results/TraitMaps/primat_Invert_sd.svg",
       width = 2, height = 1.2)
# Diet.Vertebrate
Vert <- calc_trait("Diet.Vertebrate", trait_df, comm_mat, div_metrics)
map_trait(Vert, "mean_trait")
ggsave("Results/TraitMaps/primat_Vert_mean.svg",
       width = 2, height = 1.2)
map_trait(Vert, "sd_trait")
ggsave("Results/TraitMaps/primat_Vert_sd.svg",
       width = 2, height = 1.2)
# Diet.Fruit
Fruit <- calc_trait("Diet.Fruit", trait_df, comm_mat, div_metrics)
map_trait(Fruit, "mean_trait")
ggsave("Results/TraitMaps/primat_Fruit_mean.svg",
       width = 2, height = 1.2)
map_trait(Fruit, "sd_trait")
ggsave("Results/TraitMaps/primat_Fruit_sd.svg",
       width = 2, height = 1.2)
# Diet.Leaf
Leaf <- calc_trait("Diet.Leaf", trait_df, comm_mat, div_metrics)
map_trait(Leaf, "mean_trait")
ggsave("Results/TraitMaps/primat_Leaf_mean.svg",
       width = 2, height = 1.2)
map_trait(Leaf, "sd_trait")
ggsave("Results/TraitMaps/primat_Leaf_sd.svg",
       width = 2, height = 1.2)

# Artiodactyla ----------
trait_df <- read.csv("Data/traits_artiod.csv")
comm_mat <- readRDS("Data/comm_mat_artiod.RDS")
div_metrics <- artiod_div_metrics
# logMass
logMass <- calc_trait("logMass", trait_df, comm_mat, div_metrics)
map_trait(logMass, "mean_trait")
ggsave("Results/TraitMaps/artiod_logMass_mean.svg",
       width = 2, height = 1.2)
map_trait(logMass, "sd_trait")
ggsave("Results/TraitMaps/artiod_logMass_sd.svg",
       width = 2, height = 1.2)
# Diet.Graminoids
Gram <- calc_trait("Diet.Graminoids", trait_df, comm_mat, div_metrics)
map_trait(Gram, "mean_trait")
ggsave("Results/TraitMaps/artiod_Gram_mean.svg",
       width = 2, height = 1.2)
map_trait(Gram, "sd_trait")
ggsave("Results/TraitMaps/artiod_Gram_sd.svg",
       width = 2, height = 1.2)
# Diet.Browse.Fruit
BroFru <- calc_trait("Diet.Browse.Fruit", trait_df, comm_mat, div_metrics)
map_trait(BroFru, "mean_trait")
ggsave("Results/TraitMaps/artiod_BroFru_mean.svg",
       width = 2, height = 1.2)
map_trait(BroFru, "sd_trait")
ggsave("Results/TraitMaps/artiod_BroFru_sd.svg",
       width = 2, height = 1.2)
# Diet.Meat
Meat <- calc_trait("Diet.Meat", trait_df, comm_mat, div_metrics)
map_trait(Meat, "mean_trait")
ggsave("Results/TraitMaps/artiod_Meat_mean.svg",
       width = 2, height = 1.2)
map_trait(Meat, "sd_trait")
ggsave("Results/TraitMaps/artiod_Meat_sd.svg",
       width = 2, height = 1.2)

# Carnivora ----------
trait_df <- read.csv("Data/traits_carniv.csv")
comm_mat <- readRDS("Data/comm_mat_carniv.RDS")
div_metrics <- carniv_div_metrics
# logMass
logMass <- calc_trait("logMass", trait_df, comm_mat, div_metrics)
map_trait(logMass, "mean_trait")
ggsave("Results/TraitMaps/carniv_logMass_mean.svg",
       width = 2, height = 1.2)
map_trait(logMass, "sd_trait")
ggsave("Results/TraitMaps/carniv_logMass_sd.svg",
       width = 2, height = 1.2)
# ForStrat.Arboreal
Arboreal <- calc_trait("ForStrat.Arboreal", trait_df, comm_mat, div_metrics)
map_trait(Arboreal, "mean_trait")
ggsave("Results/TraitMaps/carniv_Arboreal_mean.svg",
       width = 2, height = 1.2)
map_trait(Arboreal, "sd_trait")
ggsave("Results/TraitMaps/carniv_Arboreal_sd.svg",
       width = 2, height = 1.2)
# ForStrat.Ground
Ground <- calc_trait("ForStrat.Ground", trait_df, comm_mat, div_metrics)
map_trait(Ground, "mean_trait")
ggsave("Results/TraitMaps/carniv_Ground_mean.svg",
       width = 2, height = 1.2)
map_trait(Ground, "sd_trait")
ggsave("Results/TraitMaps/carniv_Ground_sd.svg",
       width = 2, height = 1.2)
# Diet.Invertebrate
Invert <- calc_trait("Diet.Invertebrate", trait_df, comm_mat, div_metrics)
map_trait(Invert, "mean_trait")
ggsave("Results/TraitMaps/carniv_Invert_mean.svg",
       width = 2, height = 1.2)
map_trait(Invert, "sd_trait")
ggsave("Results/TraitMaps/carniv_Invert_sd.svg",
       width = 2, height = 1.2)
# Diet.Vertebrate
Vert <- calc_trait("Diet.Vertebrate", trait_df, comm_mat, div_metrics)
map_trait(Vert, "mean_trait")
ggsave("Results/TraitMaps/carniv_Vert_mean.svg",
       width = 2, height = 1.2)
map_trait(Vert, "sd_trait")
ggsave("Results/TraitMaps/carniv_Vert_sd.svg",
       width = 2, height = 1.2)
# Diet.Plant
Plant <- calc_trait("Diet.Plant", trait_df, comm_mat, div_metrics)
map_trait(Plant, "mean_trait")
ggsave("Results/TraitMaps/carniv_Plant_mean.svg",
       width = 2, height = 1.2)
map_trait(Plant, "sd_trait")
ggsave("Results/TraitMaps/carniv_Plant_sd.svg",
       width = 2, height = 1.2)


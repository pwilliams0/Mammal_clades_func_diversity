#=========================================================================#
# Create assemblages                                                      #
# Process trait data for each order                                       #
#   (some species will need trait data manually entered before next step) #
#=========================================================================#
library(sf)
library(terra)
library(tidyverse)
library(exactextractr)

# ----- CREATE ASSEMBLAGES, PRESENT_NATURAL ----------

# Load grid cells
cells <- st_read("Data/cells.shp")

# Load all present-natural species ranges
rast_list <- list.files(path = "Data/Raw/PHYLACINE_1.2/Present_natural",
                        pattern='.tif$', all.files= T,
                        full.names= T)
ranges <- terra::rast(rast_list)

# Get list of species in each order
# RODENTIA
spp_rodent <- read.csv("Data/Raw/PHYLACINE_1.2/Trait_data.csv") %>%
  filter(Order.1.2 == "Rodentia") %>% 
  pull(Binomial.1.2)
# CHIROPTERA
spp_chirop <- read.csv("Data/Raw/PHYLACINE_1.2/Trait_data.csv") %>%
  filter(Order.1.2 == "Chiroptera") %>% 
  pull(Binomial.1.2)
# EULIPOTYPHLA
spp_eulipo <- read.csv("Data/Raw/PHYLACINE_1.2/Trait_data.csv") %>%
  filter(Order.1.2 == "Eulipotyphla") %>% 
  # Exclude Caribbean families
  filter(!Family.1.2 %in% c("Nesophontidae","Solenodontidae")) %>%
  pull(Binomial.1.2)
# PRIMATES
spp_primat <- read.csv("Data/Raw/PHYLACINE_1.2/Trait_data.csv") %>%
  filter(Order.1.2 == "Primates") %>%
  # Exclude genus Homo
  filter(Genus.1.2 != "Homo") %>%
  pull(Binomial.1.2)
# CARNIVORA
spp_carniv <- read.csv("Data/Raw/PHYLACINE_1.2/Trait_data.csv") %>%
  filter(Order.1.2 == "Carnivora") %>% 
  # Exclude pinnipeds
  filter(!Family.1.2 %in% c("Odobenidae","Otariidae","Phocidae")) %>%
  pull(Binomial.1.2)
# ARTIODACTYLA
spp_artiod <- read.csv("Data/Raw/PHYLACINE_1.2/Trait_data.csv") %>%
  filter(Order.1.2 == "Cetartiodactyla") %>% 
  # Exclude cetaceans
  filter(Terrestrial == 1) %>%
  pull(Binomial.1.2)

# Create a community matrix
# Arguments:
#   spp_list = list of species to use
#   ranges = stack of ranges (SpatRaster) to use
#   cells = shapefile of grid cells with "cell_id" column
make_comm_mat <- function(spp_list, ranges, cells) {
  # Keep ranges of species in the Order
  ranges_spp <- ranges[[spp_list]]
  # Extract names of species occurring in each cell
  # Count any with >= 1 raster cell in each grid cell (4 raster cells per grid cell)
  comm_mat <- exactextractr::exact_extract(ranges_spp, cells, 'sum')
  comm_mat[comm_mat<1] <- 0
  comm_mat[comm_mat>=1] <- 1
  colnames(comm_mat) <- names(ranges_spp)
  rownames(comm_mat) <- cells$cell_id
  # Save community matrix
  # Keep species that are present in at least one cell
  comm_mat[,which(colSums(comm_mat) > 0)]
}

# All mammals
spp_mamm <- read.csv("Data/Raw/PHYLACINE_1.2/Trait_data.csv") %>%
  mutate(Terr_Aer = Terrestrial + Aerial) %>%
  filter(Terr_Aer == 1) %>% 
  pull(Binomial.1.2)
comm_mat <- make_comm_mat(spp_mamm, ranges, cells)
saveRDS(comm_mat, "Data/comm_mat_mammals.RDS")
# Rodentia
comm_mat <- make_comm_mat(spp_rodent, ranges, cells)
saveRDS(comm_mat, "Data/comm_mat_rodent.RDS")
# Chiroptera
comm_mat <- make_comm_mat(spp_chirop, ranges, cells)
saveRDS(comm_mat, "Data/comm_mat_chirop.RDS")
# Eulipotyphla
comm_mat <- make_comm_mat(spp_eulipo, ranges, cells)
saveRDS(comm_mat, "Data/comm_mat_eulipo.RDS")
# Primates
comm_mat <- make_comm_mat(spp_primat, ranges, cells)
saveRDS(comm_mat, "Data/comm_mat_primat.RDS")
# Artiodactyla
comm_mat <- make_comm_mat(spp_artiod, ranges, cells)
saveRDS(comm_mat, "Data/comm_mat_artiod.RDS")
# Carnivora
comm_mat <- make_comm_mat(spp_carniv, ranges, cells)
saveRDS(comm_mat, "Data/comm_mat_carniv.RDS")

# ----- PROCESS TRAIT DATA ----------

# Load Phylacine & EltonTraits names to align taxonomies
elton_names <- read.csv("Data/Raw/PHYLACINE_1.2/Synonymy_table_valid_species_only.csv") %>%
  mutate(Scientific = paste(EltonTraits.1.0.Genus,
                            EltonTraits.1.0.Species, sep=" ")) %>%
  select(Binomial.1.2, Scientific)

# Load HerbiTraits, get Mass, Diet, & Habitat strata
#   Remove subspecies
HerbiTraits <- read.csv("Data/Raw/HerbiTraits_1.2.csv") %>%
  filter(Class == "Mammalia",
         Domesticated == 0) %>%
  separate(Binomial, into=c("Genus.1.2","Species.1.2")) %>%
  mutate(Mass.Herb = Mass.g) %>%
  select(Genus.1.2, Species.1.2, Mass.Herb,
         Diet.Graminoids, Diet.Browse.Fruit, Diet.Meat,
         Habitat.Terrestrial, Habitat.Arboreal)

# Load CarniTraits, get Mass, Diet, & Hunting strata
CarniTraits <- read.csv("Data/Raw/CarniTraitsv1.5.3.csv") %>%
  mutate(Mass.Carn = Mass.g,
         Diet.I.Carn = Diet.Invertebrate,
         Diet.V.Carn = Diet.Vertebrate,
         Diet.P.Carn = Diet.Plant) %>%
  select(Binomial.1.2, Mass.Carn, Diet.I.Carn, Diet.V.Carn, Diet.P.Carn,
         Hunting.terrestrial, Hunting.arboreal) %>%
  # In cases where diet categories sum to > 100 (e.g., Dusicyon_avus),
  #   set diet categories to NA
  mutate(Total.Diet = Diet.I.Carn + Diet.V.Carn + Diet.P.Carn,
         Diet.I.Carn = ifelse(Total.Diet>100, NA, Diet.I.Carn),
         Diet.V.Carn = ifelse(Total.Diet>100, NA, Diet.V.Carn),
         Diet.P.Carn = ifelse(Total.Diet>100, NA, Diet.P.Carn)) %>%
  select(-Total.Diet)

# Load EltonTraits for all mammal species
# Process trait data for all mammal species
traits_all <- read.table(file="Data/Raw/MamFuncDat.txt", sep="\t", header=TRUE,
                         stringsAsFactors = FALSE) %>%
  left_join(elton_names, "Scientific") %>%
  # Convert ForStrat variables to Arboreal/Terrestrial/Aerial
  #   Scansorial species count as both Arboreal and Terrestrial
  mutate(ForStrat.Arboreal = case_when(ForStrat.Value == "Ar" ~ 1,
                                       ForStrat.Value == "S" ~ 1,
                                       ForStrat.Value == "G" ~ 0,
                                       ForStrat.Value == "M" ~ 0,
                                       ForStrat.Value == "A" ~ 0),
         ForStrat.Ground = case_when(ForStrat.Value == "Ar" ~ 0,
                                     ForStrat.Value == "S" ~ 1,
                                     ForStrat.Value == "G" ~ 1,
                                     ForStrat.Value == "M" ~ 1,
                                     ForStrat.Value == "A" ~ 0),
         ForStrat.Aerial = ifelse(ForStrat.Value == "A", 1, 0)) %>%
  select(Binomial.1.2, Diet.Fruit, Diet.Nect, Diet.Seed, Diet.PlantO,
         ForStrat.Arboreal, ForStrat.Ground, ForStrat.Aerial) %>%
  na.omit() %>%
  # Join EltonTraits with Phylacine, keeping ForStrat and plant Diets from EltonTraits
  right_join(read.csv("Data/Raw/PHYLACINE_1.2/Trait_data.csv"),
             by = "Binomial.1.2") %>%
  select(Binomial.1.2, Order.1.2, Family.1.2, Genus.1.2,
         Species.1.2, Diet.Plant, Diet.Vertebrate, Diet.Invertebrate,
         Mass.g, ForStrat.Arboreal, ForStrat.Ground, ForStrat.Aerial,
         Diet.Fruit, Diet.Nect, Diet.Seed, Diet.PlantO) %>%
  # Join with HerbiTraits
  left_join(HerbiTraits, by = join_by(Genus.1.2, Species.1.2)) %>%
  # Join with CarniTraits
  left_join(CarniTraits, by = "Binomial.1.2") %>%
  # Use updated body mass data from HerbiTraits or CarniTraits when available
  mutate(Mass.g = case_when(Mass.Herb >= 0 ~ Mass.Herb,
                            Mass.Carn >= 0 ~ Mass.Carn,
                            .default = Mass.g),
         Mass.Source = case_when(Mass.Herb >= 0 ~ "HerbiTraits",
                                 Mass.Carn >= 0 ~ "CarniTraits",
                                 Mass.g >= 0 ~ "PHYLACINE",
                                 .default = NA)) %>%
  # Use Arboreal & Ground foraging strata data from HerbiTraits or CarniTraits when available
  mutate(ForStrat.Arboreal = case_when(Habitat.Arboreal >= 0 ~ Habitat.Arboreal,
                                       Hunting.arboreal >= 0 ~ Hunting.arboreal,
                                       .default = ForStrat.Arboreal),
         ForStrat.Ground = case_when(Habitat.Terrestrial >= 0 ~ Habitat.Terrestrial,
                                     Hunting.terrestrial >= 0 ~ Hunting.terrestrial,
                                     .default = ForStrat.Ground),
         ForStrat.Source = case_when(Habitat.Terrestrial >= 0 ~ "HerbiTraits",
                                     Hunting.terrestrial >= 0 ~ "CarniTraits",
                                     ForStrat.Ground >= 0 ~ "EltonTraits",
                                     .default = NA)) %>%
  # Use updated diet data from CarniTraits when available
  mutate(Diet.Vertebrate = ifelse(is.na(Diet.V.Carn),
                                  Diet.Vertebrate, Diet.V.Carn),
         Diet.Invertebrate = ifelse(is.na(Diet.I.Carn),
                                    Diet.Invertebrate, Diet.I.Carn),
         Diet.Plant = ifelse(is.na(Diet.P.Carn),
                             Diet.Plant, Diet.P.Carn),
         Diet.Source = case_when(Diet.V.Carn >= 0 ~ "CarniTraits",
                                 Diet.Vertebrate >= 0 ~ "PHYLACINE",
                                 .default = NA)) %>%
  # Impute ForStrat from congeners
  #   Identify ForStrat among majority of congeners
  group_by(Genus.1.2) %>%
  mutate(ForStrat.Arboreal = ifelse(is.na(ForStrat.Arboreal),
                                    round(mean(na.omit(ForStrat.Arboreal))),
                                    ForStrat.Arboreal),
         ForStrat.Ground = ifelse(is.na(ForStrat.Ground),
                                  round(mean(na.omit(ForStrat.Ground))),
                                  ForStrat.Ground),
         ForStrat.Aerial = ifelse(is.na(ForStrat.Aerial),
                                  round(mean(na.omit(ForStrat.Aerial))),
                                  ForStrat.Aerial)) %>%
  ungroup() %>%
  mutate(ForStrat.Source = case_when(!is.na(ForStrat.Source) ~ ForStrat.Source,
                                     is.na(ForStrat.Ground) ~ NA,
                                     .default = "Imputed from genus")) %>%
  # Select and order columns
  select(Binomial.1.2, Order.1.2, Family.1.2, Genus.1.2, Species.1.2,
         Mass.g, Mass.Source,
         Diet.Plant, Diet.Vertebrate, Diet.Invertebrate, Diet.Source,
         ForStrat.Arboreal, ForStrat.Ground, ForStrat.Aerial, ForStrat.Source,
         Diet.Fruit, Diet.Nect, Diet.Seed, Diet.PlantO,
         Diet.Graminoids, Diet.Browse.Fruit, Diet.Meat)

# ---------- Rodentia ---------

# Diet categories: Invertebrate, Vertebrate, Fruit, Seed, Herb
#    Nect & Fruit combined into "Fruit"
#    PlantO renamed "Herb"

# Process trait data
traits_rodent <- traits_all %>%
  # Filter for species in this order in this study
  filter(Binomial.1.2 %in% colnames(readRDS("Data/comm_mat_rodent.RDS"))) %>%
  # Calculate % of total plant diet for Fruit/Nect, Seed, & Herb
  #   If total plant diet is 0, set percentages to 0
  mutate(Total.Plant = Diet.Fruit + Diet.Nect + Diet.Seed + Diet.PlantO,
         Perc.Fruit = ifelse(Total.Plant > 0,
                                 (Diet.Fruit+Diet.Nect) / Total.Plant, 0),
         Perc.Seed = ifelse(Total.Plant > 0,
                            Diet.Seed / Total.Plant, 0),
         Perc.Herb = ifelse(Total.Plant > 0,
                            Diet.PlantO / Total.Plant, 0),
         PlantCat.Source = case_when(Diet.Plant==0 ~ "Non-Herbivore",
                                     Total.Plant > 0 ~ "EltonTraits",
                                     .default = NA)) %>%
  # Impute % plant diet categories from genus
  #   Calculate average % among species in genus that consume plants according to EltonTraits
  #   If no plant-consuming congeners, leave as NA
  #   If the species eats no plants according to EltonTraits, impute % from congeners
  #     Some species have Diet.Plant > 0 in Phylacine but no plant diet in EltonTraits,
  #     so need to impute plant categories
  #   If the species lacks EltonTraits data, impute from congeners
  group_by(Genus.1.2) %>%
  mutate(N_Plant = sum(na.omit(Total.Plant > 0)), # Number of plant-consuming species in genus
         Perc.Fruit = case_when(N_Plant == 0 ~ NA,
                                    Total.Plant == 0 ~ sum(na.omit(Perc.Fruit)) / N_Plant,
                                    is.na(Perc.Fruit) ~ sum(na.omit(Perc.Fruit)) / N_Plant,
                                    .default = Perc.Fruit),
         Perc.Seed = case_when(N_Plant == 0 ~ NA,
                               Total.Plant == 0 ~ sum(na.omit(Perc.Seed)) / N_Plant,
                               is.na(Perc.Seed) ~ sum(na.omit(Perc.Seed)) / N_Plant,
                               .default = Perc.Seed),
         Perc.Herb = case_when(N_Plant == 0 ~ NA,
                               Total.Plant == 0 ~ sum(na.omit(Perc.Herb)) / N_Plant,
                               is.na(Perc.Herb) ~ sum(na.omit(Perc.Herb)) / N_Plant,
                               .default = Perc.Herb),
         PlantCat.Source = ifelse(is.na(PlantCat.Source) & Perc.Herb >= 0,
                                  "Imputed from genus",
                                  PlantCat.Source)) %>%
  select(-N_Plant) %>%
  # Allocate Diet.Plant to plant diet categories based on percentages
  #   If Diet.Plant is 0, set diet category to 0
  mutate(Diet.Fruit = ifelse(Diet.Plant == 0, 0,
                                 round(Diet.Plant * Perc.Fruit)),
         Diet.Seed = ifelse(Diet.Plant == 0, 0,
                            round(Diet.Plant * Perc.Seed)),
         Diet.Herb = ifelse(Diet.Plant == 0, 0,
                            round(Diet.Plant * Perc.Herb))) %>%
  # Select and order columns
  select(Binomial.1.2, Order.1.2, Family.1.2, Genus.1.2, Species.1.2,
         Mass.g, Mass.Source,
         ForStrat.Arboreal, ForStrat.Ground, ForStrat.Source,
         Diet.Invertebrate, Diet.Vertebrate,
         Diet.Fruit, Diet.Seed, Diet.Herb, Diet.Source,
         PlantCat.Source, Perc.Fruit, Perc.Seed, Perc.Herb)

# Check how many species to manually add plant diet data
sum(is.na(traits_rodent$PlantCat.Source)) # 11 species
# Check how many species to manually add foraging strata data
sum(is.na(traits_rodent$ForStrat.Source)) #  8 species

# Save file for manual processing
write.csv(traits_rodent, "Data/Traits_tmp/traits_rodent_tmp.csv")

# ---------- Chiroptera ---------

# Diet categories: Invertebrate, Vertebrate, Fruit, Nect
#    Fruit, Seed, & PlantO combined into "Fruit"

# Process trait data
traits_chirop <- traits_all %>%
  # Filter for species in this order in this study
  filter(Binomial.1.2 %in% colnames(readRDS("Data/comm_mat_chirop.RDS"))) %>%
  # Calculate % of total plant diet for Fruit/Seed/PlantO & Nect
  #   If total plant diet is 0, set percentages to 0
  mutate(Total.Plant = Diet.Fruit + Diet.Nect + Diet.Seed + Diet.PlantO,
         Perc.Fruit = ifelse(Total.Plant > 0,
                             (Diet.Fruit+Diet.Seed+Diet.PlantO) / Total.Plant, 0),
         Perc.Nect = ifelse(Total.Plant > 0,
                            Diet.Nect / Total.Plant, 0),
         PlantCat.Source = case_when(Diet.Plant==0 ~ "Non-Herbivore",
                                     Total.Plant > 0 ~ "EltonTraits",
                                     .default = NA)) %>%
  # Impute % plant diet categories from genus
  #   Calculate average % among species in genus that consume plants according to EltonTraits
  #   If no plant-consuming congeners, leave as NA
  #   If the species eats no plants according to EltonTraits, impute % from congeners
  #     Some species have Diet.Plant > 0 in Phylacine but no plant diet in EltonTraits,
  #     so need to impute plant categories
  #   If the species lacks EltonTraits data, impute from congeners
  group_by(Genus.1.2) %>%
  mutate(N_Plant = sum(na.omit(Total.Plant > 0)), # Number of plant-consuming species in genus
         Perc.Fruit = case_when(N_Plant == 0 ~ NA,
                                    Total.Plant == 0 ~ sum(na.omit(Perc.Fruit)) / N_Plant,
                                    is.na(Perc.Fruit) ~ sum(na.omit(Perc.Fruit)) / N_Plant,
                                    .default = Perc.Fruit),
         Perc.Nect = case_when(N_Plant == 0 ~ NA,
                               Total.Plant == 0 ~ sum(na.omit(Perc.Nect)) / N_Plant,
                               is.na(Perc.Nect) ~ sum(na.omit(Perc.Nect)) / N_Plant,
                               .default = Perc.Nect),
         PlantCat.Source = ifelse(is.na(PlantCat.Source) & Perc.Nect >= 0,
                                  "Imputed from genus",
                                  PlantCat.Source)) %>%
  select(-N_Plant) %>%
  # For species still missing plant category data and with Diet.Plant < 10,
  #   set Diet.Plant to 0 and proportionally reallocate diet % to Vert and Inv
  # For these species which consume very little plant matter, it is not worth
  #   searching the literature to determine the types of plants that they consume
  mutate(low_plant = ifelse(Diet.Plant > 0 & Diet.Plant < 10 & is.na(Perc.Nect),
                            1, 0),
         Diet.Invertebrate = ifelse(low_plant == 1,
                                    Diet.Invertebrate + round(Diet.Plant * (Diet.Invertebrate / 100)),
                                    Diet.Invertebrate),
         Diet.Vertebrate = ifelse(low_plant == 1,
                                  Diet.Vertebrate + round(Diet.Plant * (Diet.Vertebrate / 100)),
                                  Diet.Vertebrate),
         Diet.Plant = ifelse(low_plant == 1, 0, Diet.Plant),
         PlantCat.Source = case_when(low_plant == 1 ~ "Diet.Plant < 10, reallocated to Vert and Inv",
                                     !is.na(PlantCat.Source) ~ PlantCat.Source,
                                     .default = NA)) %>%
  select(-low_plant) %>%
  # Allocate Diet.Plant to plant diet categories based on percentages
  #   If Diet.Plant is 0, set diet category to 0
  mutate(Diet.Fruit = ifelse(Diet.Plant == 0, 0,
                                 round(Diet.Plant * Perc.Fruit)),
         Diet.Nect = ifelse(Diet.Plant == 0, 0,
                            round(Diet.Plant * Perc.Nect))) %>%
  # Select and order columns
  select(Binomial.1.2, Order.1.2, Family.1.2, Genus.1.2, Species.1.2,
         Mass.g, Mass.Source,
         ForStrat.Arboreal, ForStrat.Ground, ForStrat.Aerial, ForStrat.Source,
         Diet.Invertebrate, Diet.Vertebrate,
         Diet.Fruit, Diet.Nect, Diet.Source,
         PlantCat.Source, Perc.Fruit, Perc.Nect)

# Check how many species to manually add plant diet data
sum(is.na(traits_chirop$PlantCat.Source)) # 3 species
# Check how many species to manually add foraging strata data
sum(is.na(traits_chirop$ForStrat.Source)) # 1 species

# Save file for manual processing
write.csv(traits_chirop, "Data/Traits_tmp/traits_chirop_tmp.csv")

# ---------- Eulipotyphla ---------

# Diet categories: Invertebrate, Vertebrate, Plant

# Process trait data
traits_eulipo <- traits_all %>%
  # Filter for species in this order in this study
  filter(Binomial.1.2 %in% colnames(readRDS("Data/comm_mat_eulipo.RDS"))) %>%
  select(Binomial.1.2, Order.1.2, Family.1.2, Genus.1.2, Species.1.2,
         Mass.g, Mass.Source,
         ForStrat.Arboreal, ForStrat.Ground, ForStrat.Source,
         Diet.Invertebrate, Diet.Vertebrate, Diet.Plant, Diet.Source)

# Check how many species to manually add foraging strata data
sum(is.na(traits_eulipo$ForStrat.Source)) # No species

# Save file for manual processing
write.csv(traits_eulipo, "Data/Traits_tmp/traits_eulipo_tmp.csv")

# ---------- Primates ---------

# Diet categories: Invertebrate, Vertebrate, Fruit, Leaf
#    Fruit, Nect, & Seed combined into "Fruit"
#    PlantO renamed as "Leaf"; most herbivorous non-frugivorous primates are folivores 

# Process trait data
traits_primat <- traits_all %>%
  # Filter for species in this order in this study
  filter(Binomial.1.2 %in% colnames(readRDS("Data/comm_mat_primat.RDS"))) %>%
  # Calculate % of total plant diet for Fruit (Fruit/Nect/Seed), & Leaf (PlantO)
  #   If total plant diet is 0, set percentages to 0
  mutate(Total.Plant = Diet.Fruit + Diet.Nect + Diet.Seed + Diet.PlantO,
         Perc.Fruit = ifelse(Total.Plant > 0,
                             (Diet.Fruit+Diet.Nect+Diet.Seed) / Total.Plant, 0),
         Perc.Leaf = ifelse(Total.Plant > 0,
                            Diet.PlantO / Total.Plant, 0),
         PlantCat.Source = case_when(Diet.Plant==0 ~ "Non-Herbivore",
                                     Total.Plant > 0 ~ "EltonTraits",
                                     .default = NA)) %>%
  # Impute % plant diet categories from genus
  #   Calculate average % among species in genus that consume plants according to EltonTraits
  #   If no plant-consuming congeners, leave as NA
  #   If the species eats no plants according to EltonTraits, impute % from congeners
  #     Some species have Diet.Plant > 0 in Phylacine but no plant diet in EltonTraits,
  #     so need to impute plant categories
  #   If the species lacks EltonTraits data, impute from congeners
  group_by(Genus.1.2) %>%
  mutate(N_Plant = sum(na.omit(Total.Plant > 0)), # Number of plant-consuming species in genus
         Perc.Fruit = case_when(N_Plant == 0 ~ NA,
                                    Total.Plant == 0 ~ sum(na.omit(Perc.Fruit)) / N_Plant,
                                    is.na(Perc.Fruit) ~ sum(na.omit(Perc.Fruit)) / N_Plant,
                                    .default = Perc.Fruit),
         Perc.Leaf = case_when(N_Plant == 0 ~ NA,
                               Total.Plant == 0 ~ sum(na.omit(Perc.Leaf)) / N_Plant,
                               is.na(Perc.Leaf) ~ sum(na.omit(Perc.Leaf)) / N_Plant,
                               .default = Perc.Leaf),
         PlantCat.Source = ifelse(is.na(PlantCat.Source) & Perc.Leaf >= 0,
                                  "Imputed from genus",
                                  PlantCat.Source)) %>%
  select(-N_Plant) %>%
  # Allocate Diet.Plant to plant diet categories based on percentages
  #   If Diet.Plant is 0, set diet category to 0
  mutate(Diet.Fruit = ifelse(Diet.Plant == 0, 0,
                                 round(Diet.Plant * Perc.Fruit)),
         Diet.Leaf = ifelse(Diet.Plant == 0, 0,
                            round(Diet.Plant * Perc.Leaf))) %>%
  # Select and order columns
  select(Binomial.1.2, Order.1.2, Family.1.2, Genus.1.2, Species.1.2,
         Mass.g, Mass.Source,
         ForStrat.Arboreal, ForStrat.Ground, ForStrat.Source,
         Diet.Invertebrate, Diet.Vertebrate,
         Diet.Fruit, Diet.Leaf, Diet.Source,
         PlantCat.Source, Perc.Fruit, Perc.Leaf)

# Check how many species to manually add plant diet data
sum(is.na(traits_primat$PlantCat.Source)) # 19 species
# Check how many species to manually add foraging strata data
sum(is.na(traits_primat$ForStrat.Source)) #  2 species

# Save file for manual processing
write.csv(traits_primat, "Data/Traits_tmp/traits_primat_tmp.csv")

# ---------- Artiodactyla ---------

# Plant diet categories: Graminoids, Browse.Fruit, Meat

# Process trait data
traits_artiod <- traits_all %>%
  # Filter for species in this order in this study
  filter(Binomial.1.2 %in% colnames(readRDS("Data/comm_mat_artiod.RDS"))) %>%
  # For species too small to be included in HerbiTraits,
  #   impute Diet.Graminoids, Diet.Browse.Fruit, & Diet.Meat from genus when possible
  group_by(Genus.1.2) %>%
  mutate(N_HerbiTrait = length(na.omit(Diet.Meat)), # Number of species in genus with HerbiTrait data
         Diet.Source = case_when(Diet.Meat >= 0 ~ "HerbiTraits",
                                 N_HerbiTrait > 0 ~ "Imputed from genus using HerbiTraits",
                                 .default = NA),
         Diet.Graminoids = case_when(N_HerbiTrait == 0 ~ NA,
                                     is.na(Diet.Graminoids) ~ round(mean(na.omit(Diet.Graminoids))),
                                     .default = Diet.Graminoids),
         Diet.Browse.Fruit = case_when(N_HerbiTrait == 0 ~ NA,
                                       is.na(Diet.Browse.Fruit) ~ round(mean(na.omit(Diet.Browse.Fruit))),
                                       .default = Diet.Browse.Fruit),
         Diet.Meat = case_when(N_HerbiTrait == 0 ~ NA,
                               is.na(Diet.Meat) ~ round(mean(na.omit(Diet.Meat))),
                               .default = Diet.Meat)) %>%
  ungroup() %>%
  select(-N_HerbiTrait) %>%
  # Select and order columns
  select(Binomial.1.2, Order.1.2, Family.1.2, Genus.1.2, Species.1.2,
         Mass.g, Mass.Source,
         ForStrat.Arboreal, ForStrat.Ground, ForStrat.Source,
         Diet.Graminoids, Diet.Browse.Fruit, Diet.Meat, Diet.Source)

# Check how many species to manually add diet data
sum(is.na(traits_artiod$Diet.Source)) #     18 species
# Check how many species to manually add foraging strata data
sum(is.na(traits_artiod$ForStrat.Source)) # No species

# Save file for manual processing
write.csv(traits_artiod, "Data/Traits_tmp/traits_artiod_tmp.csv")

# ---------- Carnivora ---------

# Diet categories: Invertebrate, Vertebrate, Plant

# Process trait data
traits_carniv <- traits_all %>%
  # Filter for species in this order in this study
  filter(Binomial.1.2 %in% colnames(readRDS("Data/comm_mat_carniv.RDS"))) %>%
  select(Binomial.1.2, Order.1.2, Family.1.2, Genus.1.2, Species.1.2,
         Mass.g, Mass.Source, ForStrat.Arboreal, ForStrat.Ground, ForStrat.Source,
         Diet.Invertebrate, Diet.Vertebrate, Diet.Plant, Diet.Source)

# Check how many species to manually add foraging strata data
sum(is.na(traits_carniv$ForStrat.Source)) # 1 species

# Save file for manual processing
write.csv(traits_carniv, "Data/Traits_tmp/traits_carniv_tmp.csv")

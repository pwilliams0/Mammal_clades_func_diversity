#========================================================#
# Create grid cells, create global environment variables #
#========================================================#

library(tidyverse)
library(sf)
library(terra)
library(raster)
library(exactextractr)
sf_use_s2(FALSE)

# ----- CREATE GRID CELLS ----------

# Load Phylacine raster to get spatial extent and projection
human_rast <- terra::rast("Data/Raw/PHYLACINE_1.2/Current/Homo_sapiens.tif")

# Load land polygons
land <- read_sf("Data/Raw/ne_50m_land/ne_50m_land.shp") %>%
  st_transform(st_crs(human_rast))

# Create equal-area grid
# 2 degree longitude at equator
# 193 km by 193 km at 30B0 North and 30B0 South
# Four PHYLACINE pixels per grid cell
grid <- human_rast %>% 
  st_make_grid(n = c(180,71)) %>% 
  st_as_sf() %>%
  mutate(area = st_area(.),
         cell_id = row_number())
# Keep cells >= 50% land
land_cells <- grid %>%
  st_intersection(land) %>%
  mutate(intersect_area = st_area(.),
         perc_land = as.numeric(intersect_area / area)) %>%
  dplyr::filter(perc_land >= 0.5) %>%
  pull(cell_id)
#Create cells to be used in all analyses
cells <- grid %>%
  dplyr::filter(cell_id %in% land_cells) %>%
  dplyr::mutate(cell_id = paste("X",row_number(),sep="")) %>%
  dplyr::select(-area)
st_geometry(cells) <- "geometry"

# Load Holt et al. 2013 realms
realm <- st_read("Data/Raw/newRealms/newRealms.shp") %>%
  st_wrap_dateline("WRAPDATELINE=YES") %>%
  st_transform(crs=st_crs(cells))

# Identify realm of each cell, join with cells
# When multiple realms per cell, choose larger area
realm_cells <- st_intersection(cells, realm) %>%
  mutate(polygon_area = st_area(.)) %>%
  group_by(cell_id) %>%
  arrange(desc(polygon_area)) %>%
  slice(1) %>% ungroup() %>%
  st_drop_geometry() %>%
  dplyr::select(cell_id, Realm)
cells <- cells %>%
  left_join(realm_cells, "cell_id")
# Define interior Greenland as "Palearctic"
cells$Realm[which(is.na(cells$Realm))] <- "Palearctic"

# Save grid cells
st_write(cells, "Data/cells.shp")

# ----- ENVIRONMENT, LANDMASS AREA ----------

# Manually assign each cell to a landmass
# Manually enter area of each landmass using data from Williams et al. 2024 Nature Communications
#   https://github.com/pwilliams0/Biogeography_and_global_diversity
#   (originally from Global Islands dataset)
# Count Eurasia as a single landmass
# Count North and South America as separate landmasses, divided by Isthmus of Panama
# Count Africa and Eurasia as separate landmasses, divided by Isthmus of Suez
# Save as "Landmass_area_cells.csv"

# ----- ENVIRONMENT, ELEVATION (MEAN AND RANGE) ----------

cells <- st_read("Data/cells.shp")

# GMTED2010 from USGS
elev_rast <- rast("Data/Raw/mn30_grd/w001000.adf")
# Check that raster looks correct
plot(elev_rast)

# Make cells compatible with raster
# Adjust northenmost row of cells so transformation works
st_geometry(cells)[3536:3560] <- st_geometry(cells)[3536:3560] +
  c(0,-6200) # Shifting cells 6.2 km south, 3.2% of grid length
# Set CRS, transform to equirectangular to match raster
cells <- cells %>% st_transform(crs=st_crs(elev_rast))

# Calculate elevation metrics
elev_df <- cells
# Calculate mean elevation for each grid cell
elev_df$elev_mean <- exactextractr::exact_extract(elev_rast, cells, 'mean')
# Calculate max and min elevation to then calculate elevation range
elev_df$elev_max <- exactextractr::exact_extract(elev_rast, cells, 'max')
elev_df$elev_min <- exactextractr::exact_extract(elev_rast, cells, 'min')

# Save elevation data
elev_df <- elev_df %>% 
  mutate(elev_range = elev_max - elev_min) %>%
  dplyr::select(cell_id, elev_mean, elev_range) %>%
  st_drop_geometry()
write.csv(elev_df, "Data/elev_cells.csv")

# ----- ENVIRONMENT, PRESENT CLIMATE ----------

cells <- st_read("Data/cells.shp")

# Data from worldclim.org, WorldClim 2.1
# Climate data from 1970-2000

# Load rasters of climate data, bioclimatic indicators
r01 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_1.tif")
r02 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_2.tif")
r03 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_3.tif")
r04 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_4.tif")
r05 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_5.tif")
r06 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_6.tif")
r07 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_7.tif")
r08 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_8.tif")
r09 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_9.tif")
r10 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_10.tif")
r11 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_11.tif")
r12 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_12.tif")
r13 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_13.tif")
r14 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_14.tif")
r15 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_15.tif")
r16 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_16.tif")
r17 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_17.tif")
r18 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_18.tif")
r19 <- raster("Data/Raw/wc2.1_10m_bio/wc2.1_10m_bio_19.tif")
# Create stack of all 19 bioclimatic variables, set CRS
clim_stack <- raster::stack(r01, r02, r03, r04, r05, r06, r07, r08, r09, r10,
                            r11, r12, r13, r14, r15, r16, r17, r18, r19)
crs(clim_stack) <- CRS("+proj=longlat +ellps=WGS84")

# Make cells compatible with raster
# Adjust northenmost row of cells so tranformation works
st_geometry(cells)[3536:3560] <- st_geometry(cells)[3536:3560] +
  c(0,-6200) # Shifting cells 6.2 km south, 3.2% of grid length
# Set CRS, transform to equirectangular to match raster
cells <- cells %>% st_transform(crs=st_crs(clim_stack))

# Calculate mean climate value per polygon, repeat for each raster layer
mean_r01_val <- exactextractr::exact_extract(clim_stack[[1]], cells, 'mean')
mean_r02_val <- exactextractr::exact_extract(clim_stack[[2]], cells, 'mean')
mean_r03_val <- exactextractr::exact_extract(clim_stack[[3]], cells, 'mean')
mean_r04_val <- exactextractr::exact_extract(clim_stack[[4]], cells, 'mean')
mean_r05_val <- exactextractr::exact_extract(clim_stack[[5]], cells, 'mean')
mean_r06_val <- exactextractr::exact_extract(clim_stack[[6]], cells, 'mean')
mean_r07_val <- exactextractr::exact_extract(clim_stack[[7]], cells, 'mean')
mean_r08_val <- exactextractr::exact_extract(clim_stack[[8]], cells, 'mean')
mean_r09_val <- exactextractr::exact_extract(clim_stack[[9]], cells, 'mean')
mean_r10_val <- exactextractr::exact_extract(clim_stack[[10]], cells, 'mean')
mean_r11_val <- exactextractr::exact_extract(clim_stack[[11]], cells, 'mean')
mean_r12_val <- exactextractr::exact_extract(clim_stack[[12]], cells, 'mean')
mean_r13_val <- exactextractr::exact_extract(clim_stack[[13]], cells, 'mean')
mean_r14_val <- exactextractr::exact_extract(clim_stack[[14]], cells, 'mean')
mean_r15_val <- exactextractr::exact_extract(clim_stack[[15]], cells, 'mean')
mean_r16_val <- exactextractr::exact_extract(clim_stack[[16]], cells, 'mean')
mean_r17_val <- exactextractr::exact_extract(clim_stack[[17]], cells, 'mean')
mean_r18_val <- exactextractr::exact_extract(clim_stack[[18]], cells, 'mean')
mean_r19_val <- exactextractr::exact_extract(clim_stack[[19]], cells, 'mean')

# Save as dataframe
#   sqrt transform precipitation variables
bioclim_Present_cells <- as.data.frame(cells) %>%
  dplyr::select(cell_id) %>%
  mutate(BioClim01 = mean_r01_val,
         BioClim02 = mean_r02_val,
         BioClim03 = mean_r03_val,
         BioClim04 = mean_r04_val,
         BioClim05 = mean_r05_val,
         BioClim06 = mean_r06_val,
         BioClim07 = mean_r07_val,
         BioClim08 = mean_r08_val,
         BioClim09 = mean_r09_val,
         BioClim10 = mean_r10_val,
         BioClim11 = mean_r11_val,
         BioClim12 = sqrt(mean_r12_val),
         BioClim13 = sqrt(mean_r13_val),
         BioClim14 = sqrt(mean_r14_val),
         BioClim15 = mean_r15_val,
         BioClim16 = sqrt(mean_r16_val),
         BioClim17 = sqrt(mean_r17_val),
         BioClim18 = sqrt(mean_r18_val),
         BioClim19 = sqrt(mean_r19_val))

# Run PCA on all cells
clim_pca <- prcomp(bioclim_Present_cells[,c(2:20)],
                   center = TRUE,scale. = TRUE)
summary(clim_pca)

# Keep first 4 axes (explains >90% variance)
clim_Present_cells <- bioclim_Present_cells %>%
  mutate(clim_pca_1 = clim_pca[["x"]][,1],
         clim_pca_2 = clim_pca[["x"]][,2],
         clim_pca_3 = clim_pca[["x"]][,3],
         clim_pca_4 = clim_pca[["x"]][,4]) %>%
  dplyr::select(cell_id, clim_pca_1, clim_pca_2,
                clim_pca_3, clim_pca_4)

# Save climate data
write.csv(clim_Present_cells, "Data/climate_cells.csv")

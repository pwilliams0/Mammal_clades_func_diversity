#==================================================#
# Run variance partitioning analysis, save results #
# Create maps of difference in residuals           #
#==================================================#

library(tidyverse)
library(sf)
library(terra)
library(colorspace)
library(rnaturalearth)

# ----- LOAD DATA ----------

cells <- st_read("Data/cells.shp") %>%
  # Landmass area
  left_join(read.csv("Data/Landmass_area_cells.csv"), "cell_id") %>%
  # Elevation (mean and range)
  left_join(read.csv("Data/elev_cells.csv"), "cell_id") %>%
  dplyr::select(-X) %>%
  # Climate
  left_join(read.csv("Data/climate_cells.csv"), "cell_id") %>%
  dplyr::select(-X)

# Load world outline
world <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_transform(crs=crs(cells))

# ----- RUN ANALYSIS, PLOT DIFFERENCE IN RESIDUALS ----------

# Run variance partitioning analysis & plot difference in residuals
#   Plot diff in residuals to see where phylobetadiversity improves model fit
# Arguments:
#   order = 6-letter taxonomic order abbreviation
#      ("rodent", "chirop", "eulipo", "primat", "carniv", or "artiod")
#   metric = diversity metric ("SR", "FRic", "SES.FRic", or "FDis")
#   metric_file = dataframe containing diversity metric and cell_id
# Note: "cells" and "world" files must be preloaded
analysis <- function(order, metric, metric_file) {
  
  # Load data
  data <- cells %>%
    left_join(read.csv(paste("Data/pb_NMDS_3_",order,".csv",sep="")),
              "cell_id") %>%
    dplyr::select(-X) %>%
    right_join(metric_file, "cell_id") %>%
    dplyr::select(-X) %>%
    dplyr::filter(!is.na(get(metric)))
  
  # Run model, all covariates
  model_all <- lm(get(metric) ~ log(Area_km2) +
                    elev_mean + I(elev_mean^2) +
                    elev_range +
                    clim_pca_1 + I(clim_pca_1^2) +
                    clim_pca_2 + I(clim_pca_2^2) +
                    clim_pca_3 + I(clim_pca_3^2) +
                    clim_pca_4 + I(clim_pca_4^2) +
                    MDS1 + I(MDS1^2) + MDS2 +
                    I(MDS2^2) + MDS3 + I(MDS3^2) +
                    MDS1:MDS3 + MDS1:MDS2 + MDS2:MDS3,
                  data=data,
                  na.action = na.fail)
  # Run model, environment
  model_env <- lm(get(metric) ~ log(Area_km2) +
                    elev_mean + I(elev_mean^2) +
                    elev_range +
                    clim_pca_1 + I(clim_pca_1^2) +
                    clim_pca_2 + I(clim_pca_2^2) +
                    clim_pca_3 + I(clim_pca_3^2) +
                    clim_pca_4 + I(clim_pca_4^2),
                  data=data,
                  na.action = na.fail)
  # Run model, phylobetadiversity
  model_pb <- lm(get(metric) ~ MDS1 + I(MDS1^2) + MDS2 +
                   I(MDS2^2) + MDS3 + I(MDS3^2) +
                   MDS1:MDS3 + MDS1:MDS2 + MDS2:MDS3,
                 data=data,
                 na.action = na.fail)
  # Save R2 of each model
  R2_all <- summary(model_all)$adj.r.squared
  R2_env <- summary(model_env)$adj.r.squared
  R2_pb <- summary(model_pb)$adj.r.squared
  # Phylobetadiversity only
  R2_pb_only <- R2_all - R2_env
  # Environment only
  R2_env_only <- R2_all - R2_pb
  # Shared
  R2_shared <- R2_all - R2_pb_only - R2_env_only
  
  # Save variance partitioning results
  results <- data.frame(matrix(ncol = 5, nrow = 3))
  colnames(results) <- c("Order","Metric","varpart","R2")
  results$Order <- order
  results$Metric <- metric
  results$varpart <- factor(c("Shared","Environment only",
                              "Phylobetadiversity only"),
                            levels=c("Shared","Environment only",
                                     "Phylobetadiversity only"))
  results$R2 <- c(R2_shared, R2_env_only, R2_pb_only)
  write.csv(results,
            paste("Data/Output/",metric,"_",order,".csv", sep=""))
 
  # Map differences in residuals with vs. without phylobetadiversity
  
  # Calculate differences in residuals
  res_df <- data.frame(
    res = abs(resid(model_env)) - abs(resid(model_all)),
    cell_id = data$cell_id)
  res_map <- cells %>% right_join(res_df, by="cell_id")
  
  # Create figure
  ggplot() +
    geom_sf(data = world, colour="grey85", linewidth=.1) +
    geom_sf(data = res_map, aes(color = res, fill = res),
            linewidth=.1) +
    scale_fill_continuous_diverging("Blue-Red 3") +
    scale_color_continuous_diverging("Blue-Red 3") +
    theme_void() +
    theme(legend.title = element_blank(),
          legend.text = element_text(size=7),
          legend.position = "bottom") +
    theme(legend.key.height = unit(.12, "cm"),
          legend.key.width = unit(.8,"cm"))
  
}

# ---------- Rodentia ----------

data <- read.csv("Data/rodent_div_metrics.csv")

# Species Richness
analysis("rodent", "SRic", data)
ggsave("Results/diffres_SRic_rodent.svg", width = 2, height = 1.2)

# SRic-corrected Functional Richness (SES.FRic)
analysis("rodent", "SES.FRic", data)
ggsave("Results/diffres_SES.FRic_rodent.svg", width = 2, height = 1.2)

# Functional Dispersion (FDis)
analysis("rodent", "FDis", data)
ggsave("Results/diffres_FDis_rodent.svg", width = 2, height = 1.2)

# ---------- Chiroptera ----------

data <- read.csv("Data/chirop_div_metrics.csv")

# Species Richness
analysis("chirop", "SRic", data)
ggsave("Results/diffres_SRic_chirop.svg", width = 2, height = 1.2)

# SRic-corrected Functional Richness (SES.FRic)
analysis("chirop", "SES.FRic", data)
ggsave("Results/diffres_SES.FRic_chirop.svg", width = 2, height = 1.2)

# Functional Dispersion (FDis)
analysis("chirop", "FDis", data)
ggsave("Results/diffres_FDis_chirop.svg", width = 2, height = 1.2)

# ---------- Eulipotyphla ----------

data <- read.csv("Data/eulipo_div_metrics.csv")

# Species Richness
analysis("eulipo", "SRic", data)
ggsave("Results/diffres_SRic_eulipo.svg", width = 2, height = 1.2)

# SRic-corrected Functional Richness (SES.FRic)
analysis("eulipo", "SES.FRic", data)
ggsave("Results/diffres_SES.FRic_eulipo.svg", width = 2, height = 1.2)

# Functional Dispersion (FDis)
analysis("eulipo", "FDis", data)
ggsave("Results/diffres_FDis_eulipo.svg", width = 2, height = 1.2)

# ---------- Primates ----------

data <- read.csv("Data/primat_div_metrics.csv")

# Species Richness
analysis("primat", "SRic", data)
ggsave("Results/diffres_SRic_primat.svg", width = 2, height = 1.2)

# SRic-corrected Functional Richness (SES.FRic)
analysis("primat", "SES.FRic", data)
ggsave("Results/diffres_SES.FRic_primat.svg", width = 2, height = 1.2)

# Functional Dispersion (FDis)
analysis("primat", "FDis", data)
ggsave("Results/diffres_FDis_primat.svg", width = 2, height = 1.2)

# ---------- Artiodactyla ----------

data <- read.csv("Data/artiod_div_metrics.csv")

# Species Richness
analysis("artiod", "SRic", data)
ggsave("Results/diffres_SRic_artiod.svg", width = 2, height = 1.2)

# SRic-corrected Functional Richness (SES.FRic)
analysis("artiod", "SES.FRic", data)
ggsave("Results/diffres_SES.FRic_artiod.svg", width = 2, height = 1.2)

# Functional Dispersion (FDis)
analysis("artiod", "FDis", data)
ggsave("Results/diffres_FDis_artiod.svg", width = 2, height = 1.2)

# ---------- Carnivora ----------

data <- read.csv("Data/carniv_div_metrics.csv")

# Species Richness
analysis("carniv", "SRic", data)
ggsave("Results/diffres_SRic_carniv.svg", width = 2, height = 1.2)

# SRic-corrected Functional Richness (SES.FRic)
analysis("carniv", "SES.FRic", data)
ggsave("Results/diffres_SES.FRic_carniv.svg", width = 2, height = 1.2)

# Functional Dispersion (FDis)
analysis("carniv", "FDis", data)
ggsave("Results/diffres_FDis_carniv.svg", width = 2, height = 1.2)

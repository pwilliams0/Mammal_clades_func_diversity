#=============================#
# Run phylobetadiversity NMDS #
#=============================#

library(tidyverse)
library(vegan)

# Load phylobetadiversity distance matrix
#   phylo.beta.sim is the turnover distance matrix
pb_dist_rodent <- readRDS("Data/pb_rodent.rds")$phylo.beta.sim
pb_dist_chirop <- readRDS("Data/pb_chirop.rds")$phylo.beta.sim
pb_dist_eulipo <- readRDS("Data/pb_eulipo.rds")$phylo.beta.sim
pb_dist_primat <- readRDS("Data/pb_primat.rds")$phylo.beta.sim
pb_dist_carniv <- readRDS("Data/pb_carniv.rds")$phylo.beta.sim
pb_dist_artiod <- readRDS("Data/pb_artiod.rds")$phylo.beta.sim

# Run phylobetadiversity NMDS
# Arguments:
#   pb_dist = phylobetadiversity distance matrix
#   n_axes = number of NMDS axes
#   trymax = max number of randmom starts when running NMDS
pb_nmds <- function(pb_dist, n_axes, trymax) {
  # Run NMDS
  set.seed(1234)
  phylo_mds <- vegan::metaMDS(pb_dist, k=n_axes, trymax=trymax)
  # Check stress
  print(phylo_mds$stress)
  # Save NMDS axes
  pb_NMDS <- as.data.frame(phylo_mds$points) %>%
    rownames_to_column("cell_id")
  return(pb_NMDS)
}

# ----- RUN NMDS, 2 axes ----------

# Rodentia     (stress = 0.2224307)
pb_NMDS_2 <- pb_nmds(pb_dist_rodent, n_axes=2, trymax=20)
write.csv(pb_NMDS_2, "Data/pb_NMDS_2_rodent.csv")
# Chiroptera   (stress = 0.2011103)
pb_NMDS_2 <- pb_nmds(pb_dist_chirop, n_axes=2, trymax=20)
write.csv(pb_NMDS_2, "Data/pb_NMDS_2_chirop.csv")
# Eulipotyphla (stress = 0.2169726)
pb_NMDS_2 <- pb_nmds(pb_dist_eulipo, n_axes=2, trymax=20)
write.csv(pb_NMDS_2, "Data/pb_NMDS_2_eulipo.csv")
# Primates     (stress = 0.1196383)
pb_NMDS_2 <- pb_nmds(pb_dist_primat, n_axes=2, trymax=20)
write.csv(pb_NMDS_2, "Data/pb_NMDS_2_primat.csv")
# Carnivora    (stress = 0.1461788)
pb_NMDS_2 <- pb_nmds(pb_dist_carniv, n_axes=2, trymax=20)
write.csv(pb_NMDS_2, "Data/pb_NMDS_2_carniv.csv")
# Artiodactyla (stress = 0.2087867)
pb_NMDS_2 <- pb_nmds(pb_dist_artiod, n_axes=2, trymax=20)
write.csv(pb_NMDS_2, "Data/pb_NMDS_2_artiod.csv")

# ----- RUN NMDS, 3 axes ----------

# Rodentia     (stress = 0.1635821)
pb_NMDS_3 <- pb_nmds(pb_dist_rodent, n_axes=3, trymax=20)
write.csv(pb_NMDS_3, "Data/pb_NMDS_3_rodent.csv")
# Chiroptera   (stress = 0.1583092)
pb_NMDS_3 <- pb_nmds(pb_dist_chirop, n_axes=3, trymax=20)
write.csv(pb_NMDS_3, "Data/pb_NMDS_3_chirop.csv")
# Eulipotyphla (stress = 0.1608137)
pb_NMDS_3 <- pb_nmds(pb_dist_eulipo, n_axes=3, trymax=20)
write.csv(pb_NMDS_3, "Data/pb_NMDS_3_eulipo.csv")
# Primates     (stress = 0.1082122)
pb_NMDS_3 <- pb_nmds(pb_dist_primat, n_axes=3, trymax=20)
write.csv(pb_NMDS_3, "Data/pb_NMDS_3_primat.csv")
# Carnivora    (stress = 0.1158708)
pb_NMDS_3 <- pb_nmds(pb_dist_carniv, n_axes=3, trymax=20)
write.csv(pb_NMDS_3, "Data/pb_NMDS_3_carniv.csv")
# Artiodactyla (stress = 0.1507749)
pb_NMDS_3 <- pb_nmds(pb_dist_artiod, n_axes=3, trymax=20)
write.csv(pb_NMDS_3, "Data/pb_NMDS_3_artiod.csv")


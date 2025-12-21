#==============================#
# Calculate phylobetadiversity #
#==============================#

library(picante)
library(ape)
library(betapart)

# Load community matrices, Present_natural
rodent_mat <- readRDS("Data/comm_mat_rodent.RDS")
chirop_mat <- readRDS("Data/comm_mat_chirop.RDS")
eulipo_mat <- readRDS("Data/comm_mat_eulipo.RDS")
primat_mat <- readRDS("Data/comm_mat_primat.RDS")
carniv_mat <- readRDS("Data/comm_mat_carniv.RDS")
artiod_mat <- readRDS("Data/comm_mat_artiod.RDS")

# Load phylogeny
phylo <- read.tree("Data/mamm_phylo.tre")

# Calculate phylobetadiveristy
# Arguments:
#   comm_mat = community matrix
#   phylo_all = phylogeny
calc_pb <- function(comm_mat, phylo) {
  print(sys.call())
  # Prune community matrix for cells containing at least one species
  comm_mat_subset <- comm_mat[rowSums(comm_mat) > 0,]
  # Prune phylogeny for species in community matrix
  phylo_subset <- ape::multi2di(
    picante::prune.sample(comm_mat_subset,phylo),
    random=TRUE)
  print("Done loading data")
  # Run phylogenetic beta diversity (Sorenson)
  pb <- betapart::phylo.beta.pair(
    comm_mat_subset, phylo_subset,
    index.family="sorensen")
  print("Done calculating phylobetadiversity")
  return(pb)
}

# Rodentia    (375.4 GB, 10 hr 42 min)
pb <- calc_pb(rodent_mat, phylo)
saveRDS(pb, "Data/pb_rodent.rds")
# Chiroptera  (153.6 GB,  5 hr  9 min)
pb <- calc_pb(chirop_mat, phylo)
saveRDS(pb, "Data/pb_chirop.rds")
# Eulipotyphla (44.1 GB,  1 hr 52 min)
pb <- calc_pb(eulipo_mat, phylo)
saveRDS(pb, "Data/pb_eulipo.rds")
# Primates      (8.5 GB,  0 hr 21 min)
pb <- calc_pb(primat_mat, phylo)
saveRDS(pb, "Data/pb_primat.rds")
# Carnivora    (36.1 GB,  2 hr  4 min) 
pb <- calc_pb(carniv_mat, phylo)
saveRDS(pb, "Data/pb_carniv.rds")
# Artiodactyla (41.2 GB,  1 hr 51 min)
pb <- calc_pb(artiod_mat, phylo)
saveRDS(pb, "Data/pb_artiod.rds")

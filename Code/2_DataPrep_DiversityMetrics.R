#==================================================#
# Calculate species richness (SRic),               #
#   functional richness (FRic),                    #
#   SRic-corrected functional richness (SES.FRic), #
#   & functional dispersion (FDis)                 #
#==================================================#
library(tidyverse)
library(mFD)
library(picante)
library(matrixStats)

# ----- CALCULATE FUNCTIONAL PCoA ----------

# Create PCoA using functional traits
# Arguments:
#   spp_traits = species traits, with species names as row names
#   traits_cat = trait categories
make_func_pcoa <- function(spp_traits, traits_cat) {
  # Create distance matrix of species using functional traits
  spp_dist <- mFD::funct.dist(
    sp_tr         = spp_traits,
    tr_cat        = traits_cat,
    metric        = "gower",
    scale_euclid  = "scale_center", # Standardize variables
    ordinal_var   = "classic",
    weight_type   = "equal",
    stop_if_NA    = TRUE)
  # Calculate PCoA using distance matrix
  fspaces_quality <- mFD::quality.fspaces(
    sp_dist             = spp_dist,
    maxdim_pcoa         = 5,
    deviation_weighting = c("absolute", "squared"),
    fdendro             = NULL)
  # Save functional space
  fspaces_quality
}

# ---------- Rodentia ----------

# Load trait data
traits_rodent <- read.csv("Data/traits_rodent.csv") %>%
  column_to_rownames("Binomial.1.2") %>%
  mutate(log.Mass = log(Mass.g)) %>%
  select(log.Mass, ForStrat.Arboreal, ForStrat.Ground,
         Diet.Invertebrate, Diet.Vertebrate,
         Diet.Fruit, Diet.Seed, Diet.Herb)

# Load trait categories
traits_cat <- data.frame(
  trait_name = c("log.Mass", "ForStrat.Arboreal", "ForStrat.Ground",
                 "Diet.Invertebrate", "Diet.Vertebrate",
                 "Diet.Fruit", "Diet.Seed", "Diet.Herb"),
  trait_type = c("Q", "F", "F", "F", "F", "F", "F", "F"),
  fuzzy_name = c(NA, "ForStrat", "ForStrat",
                 "Diet", "Diet", "Diet", "Diet", "Diet")) 

# Caclulate PCoA
func_space <- make_func_pcoa(traits_rodent, traits_cat)
# Check quality of PCoA
func_space$quality_fspaces
# Lowest mad & rmsd score is 4 axes, use 4 axes
func_PCoA <- func_space$details_fspaces$sp_pc_coord[,1:4]
saveRDS(func_PCoA, "Data/rodent_func_PCoA.RDS")

# # Check which traits are associated with PCoA functioal axes
# func_axes <- mFD::traits.faxes.cor(
#   sp_tr = traits_rodent, sp_faxes_coord = func_PCoA, plot = TRUE)
# View(func_axes$tr_faxes_stat)
# func_axes$"tr_faxes_plot"
# # Check functional space, pairs of functional axes
# big_plot <- funct.space.plot(func_PCoA)
# big_plot$patchwork

# ---------- Chiroptera ----------

# Load trait data
traits_chirop <- read.csv("Data/traits_chirop.csv") %>%
  column_to_rownames("Binomial.1.2") %>%
  mutate(log.Mass = log(Mass.g)) %>%
  select(log.Mass, ForStrat.Arboreal, ForStrat.Ground, ForStrat.Aerial,
         Diet.Invertebrate, Diet.Vertebrate, Diet.Fruit, Diet.Nect)

# Load trait categories
traits_cat <- data.frame(
  trait_name = c("log.Mass", "ForStrat.Arboreal", "ForStrat.Ground", "ForStrat.Aerial",
                 "Diet.Invertebrate", "Diet.Vertebrate", "Diet.Fruit", "Diet.Nect"),
  trait_type = c("Q", "F", "F", "F", "F", "F", "F", "F"),
  fuzzy_name = c(NA, "ForStrat", "ForStrat", "ForStrat",
                 "Diet", "Diet", "Diet", "Diet")) 

# Caclulate PCoA
func_space <- make_func_pcoa(traits_chirop, traits_cat)
# Check quality of PCoA
func_space$quality_fspaces
# Lowest mad & rmsd score is 4 axes, use 4 axes
func_PCoA <- func_space$details_fspaces$sp_pc_coord[,1:4]
saveRDS(func_PCoA, "Data/chirop_func_PCoA.RDS")

# # Check which traits are associated with PCoA functioal axes
# func_axes <- mFD::traits.faxes.cor(
#   sp_tr = traits_chirop, sp_faxes_coord = func_PCoA, plot = TRUE)
# View(func_axes$tr_faxes_stat)
# func_axes$"tr_faxes_plot"
# # Check functional space, pairs of functional axes
# big_plot <- funct.space.plot(func_PCoA)
# big_plot$patchwork

# ---------- Eulipotyphla ----------

# Load trait data
traits_eulipo <- read.csv("Data/traits_eulipo.csv") %>%
  column_to_rownames("Binomial.1.2") %>%
  mutate(log.Mass = log(Mass.g)) %>%
  select(log.Mass, ForStrat.Arboreal, ForStrat.Ground,
         Diet.Invertebrate, Diet.Vertebrate, Diet.Plant)

# Load trait categories
traits_cat <- data.frame(
  trait_name = c("log.Mass", "ForStrat.Arboreal", "ForStrat.Ground",
                 "Diet.Invertebrate", "Diet.Vertebrate", "Diet.Plant"),
  trait_type = c("Q", "F", "F", "F", "F", "F"),
  fuzzy_name = c(NA, "ForStrat", "ForStrat", "Diet", "Diet", "Diet")) 

# Caclulate PCoA
func_space <- make_func_pcoa(traits_eulipo, traits_cat)
# Check quality of PCoA
func_space$quality_fspaces
# Lowest nad scire is 2 axes but lowest rmsd score is 3 axes, use 3 axes
func_PCoA <- func_space$details_fspaces$sp_pc_coord[,1:3]
saveRDS(func_PCoA, "Data/eulipo_func_PCoA.RDS")

# # Check which traits are associated with PCoA functioal axes
# func_axes <- mFD::traits.faxes.cor(
#   sp_tr = traits_eulipo, sp_faxes_coord = func_PCoA, plot = TRUE)
# View(func_axes$tr_faxes_stat)
# func_axes$"tr_faxes_plot"
# # Check functional space, pairs of functional axes
# big_plot <- funct.space.plot(func_PCoA)
# big_plot$patchwork

# ---------- Primates ----------

# Load trait data
traits_primat <- read.csv("Data/traits_primat.csv") %>%
  column_to_rownames("Binomial.1.2") %>%
  mutate(log.Mass = log(Mass.g)) %>%
  select(log.Mass, ForStrat.Arboreal, ForStrat.Ground,
         Diet.Invertebrate, Diet.Vertebrate, Diet.Fruit, Diet.Leaf)

# Load trait categories
traits_cat <- data.frame(
  trait_name = c("log.Mass", "ForStrat.Arboreal", "ForStrat.Ground",
                 "Diet.Invertebrate", "Diet.Vertebrate", "Diet.Fruit", "Diet.Leaf"),
  trait_type = c("Q", "F", "F", "F", "F", "F", "F"),
  fuzzy_name = c(NA, "ForStrat", "ForStrat", "Diet", "Diet", "Diet", "Diet")) 

# Caclulate PCoA
func_space <- make_func_pcoa(traits_primat, traits_cat)
# Check quality of PCoA
func_space$quality_fspaces
# Lowest mad & rmsd score is 4 axes, use 4 axes
func_PCoA <- func_space$details_fspaces$sp_pc_coord[,1:4]
saveRDS(func_PCoA, "Data/primat_func_PCoA.RDS")

# # Check which traits are associated with PCoA functioal axes
# func_axes <- mFD::traits.faxes.cor(
#   sp_tr = traits_primat, sp_faxes_coord = func_PCoA, plot = TRUE)
# View(func_axes$tr_faxes_stat)
# func_axes$"tr_faxes_plot"
# # Check functional space, pairs of functional axes
# big_plot <- funct.space.plot(func_PCoA)
# big_plot$patchwork

# ---------- Artiodactyla ----------

# Load trait data
traits_artiod <- read.csv("Data/traits_artiod.csv") %>%
  column_to_rownames("Binomial.1.2") %>%
  mutate(log.Mass = log(Mass.g)) %>%
  select(log.Mass, Diet.Graminoids, Diet.Browse.Fruit, Diet.Meat)

# Load trait categories
traits_cat <- data.frame(
  trait_name = c("log.Mass", "Diet.Graminoids", "Diet.Browse.Fruit", "Diet.Meat"),
  trait_type = c("Q", "F", "F", "F"),
  fuzzy_name = c(NA, "Diet", "Diet", "Diet")) 

# Caclulate PCoA
func_space <- make_func_pcoa(traits_artiod, traits_cat)
# Check quality of PCoA
func_space$quality_fspaces
# Lowest mad & rmsd score is 3 axes, use 3 axes
func_PCoA <- func_space$details_fspaces$sp_pc_coord[,1:3]
saveRDS(func_PCoA, "Data/artiod_func_PCoA.RDS")

# # Check which traits are associated with PCoA functioal axes
# func_axes <- mFD::traits.faxes.cor(
#   sp_tr = traits_artiod, sp_faxes_coord = func_PCoA, plot = TRUE)
# View(func_axes$tr_faxes_stat)
# func_axes$"tr_faxes_plot"
# # Check functional space, pairs of functional axes
# big_plot <- funct.space.plot(func_PCoA)
# big_plot$patchwork
 
# ---------- Carnivora ----------

# Load trait data
traits_carniv <- read.csv("Data/traits_carniv.csv") %>%
  column_to_rownames("Binomial.1.2") %>%
  mutate(log.Mass = log(Mass.g)) %>%
  select(log.Mass, ForStrat.Arboreal, ForStrat.Ground,
         Diet.Invertebrate, Diet.Vertebrate, Diet.Plant)

# Load trait categories
traits_cat <- data.frame(
  trait_name = c("log.Mass", "ForStrat.Arboreal", "ForStrat.Ground",
                 "Diet.Invertebrate", "Diet.Vertebrate", "Diet.Plant"),
  trait_type = c("Q", "F", "F", "F", "F", "F"),
  fuzzy_name = c(NA, "ForStrat", "ForStrat", "Diet", "Diet", "Diet")) 

# Caclulate PCoA
func_space <- make_func_pcoa(traits_carniv, traits_cat)
# Check quality of PCoA
func_space$quality_fspaces
# Lowest mad & rmsd score is 3 axes, use 3 axes
func_PCoA <- func_space$details_fspaces$sp_pc_coord[,1:3]
saveRDS(func_PCoA, "Data/carniv_func_PCoA.RDS")

# # Check which traits are associated with PCoA functioal axes
# func_axes <- mFD::traits.faxes.cor(
#   sp_tr = traits_carniv, sp_faxes_coord = func_PCoA, plot = TRUE)
# View(func_axes$tr_faxes_stat)
# func_axes$"tr_faxes_plot"
# # Check functional space, pairs of functional axes
# big_plot <- funct.space.plot(func_PCoA)
# big_plot$patchwork

# ----- CALCULATE DIVERSITY METRICS ----------

# Calculate diversity metrics:
#   species richness, functional richness,
#   SES functional richness, & functional dispersion
# Arguments:
#   func_PCoA = functional PCoA
#   comm_mat = community matrix
#   iterations = number of iterations to run for SES functional richness
div_metrics <- function(comm_mat, func_PCoA, iterations) {
  
  # Calculate species richness
  SRic <- comm_mat %>%
    mutate(SRic = rowSums(.)) %>%
    select(SRic) %>%
    filter(SRic > 0) %>%
    rownames_to_column("cell_id")
  
  # Format community matrix for functional diversity calculation
  comm_matrix <- comm_mat[,rownames(func_PCoA)] %>%
    mutate(SR = rowSums(.)) %>%
    filter(SR > ncol(func_PCoA)) %>%
    dplyr::select(-SR) %>%
    as.matrix()
  
  # Calculate functional richness & functional dispersion
  FD_results <- mFD::alpha.fd.multidim(
    sp_faxes_coord = func_PCoA,
    asb_sp_w = comm_matrix,
    ind_vect = c("fric", "fdis"),
    scaling = TRUE, check_input = TRUE, verbose = TRUE,
    details_returned = TRUE)$functional_diversity_indices
  # If FRic is NA, replace with 0
  #   This can happen if there are fewer functionally distinct species than axes
  #   due to functionally identical species
  FD_results <- FD_results %>% replace(is.na(.), 0)
  
  # Calculate SES functional richness using null model
  FD_null <- matrix(nrow=nrow(FD_results), ncol=iterations)
  rownames(FD_null) <- rownames(FD_results)
  for(i in 1:iterations){
    # Create a randomized matrix controlling for species richness
    null_matrix <- randomizeMatrix(comm_matrix, null.model = "richness") 
    # Calculate FRic for randomized matrix
    temp <- mFD::alpha.fd.multidim(
      sp_faxes_coord = func_PCoA,
      asb_sp_w = null_matrix,
      ind_vect = c("fric"),
      scaling = TRUE, check_input = TRUE, verbose = FALSE,
      details_returned = TRUE)$functional_diversity_indices
    # Save FRic for this iteration
    FD_null[,i] <- temp$fric
    print(paste("Iteration",i,"of",iterations,"completed", sep=" "))
  }
  # If FRic is NA, replace with 0
  FD_null <- FD_null %>% replace(is.na(.), 0)
  
  # Save file
  div_metrics <- data.frame("cell_id" = rownames(FD_results),
                         "FRic" = FD_results$fric,
                         "FDis" = FD_results$fdis,
                         "mean_null_FRic" = rowMeans(FD_null, na.rm=TRUE),
                         "sd_null_FRic" = rowSds(FD_null, nam.rm=TRUE)) %>%
    mutate(SES.FRic = (FRic - mean_null_FRic)/sd_null_FRic) %>%
    right_join(SRic, by="cell_id") %>%
    select(cell_id, SRic, FRic, SES.FRic, FDis)
  
  div_metrics
}

# Rodentia
data <- div_metrics(comm_mat = readRDS("Data/comm_mat_rodent.RDS"),
                    func_PCoA = readRDS("Data/rodent_func_PCoA.RDS"),
                    iterations = 100)
write.csv(data, "Data/rodent_div_metrics.csv")

# Chiroptera
data <- div_metrics(comm_mat = readRDS("Data/comm_mat_chirop.RDS"),
                    func_PCoA = readRDS("Data/chirop_func_PCoA.RDS"),
                    iterations = 100)
write.csv(data, "Data/chirop_div_metrics.csv")

# Eulipotyphla
data <- div_metrics(comm_mat = readRDS("Data/comm_mat_eulipo.RDS"),
                    func_PCoA = readRDS("Data/eulipo_func_PCoA.RDS"),
                    iterations = 100)
write.csv(data, "Data/eulipo_div_metrics.csv")

# Primates
data <- div_metrics(comm_mat = readRDS("Data/comm_mat_primat.RDS"),
                    func_PCoA = readRDS("Data/primat_func_PCoA.RDS"),
                    iterations = 100)
write.csv(data, "Data/primat_div_metrics.csv")

# Artiodactyla
data <- div_metrics(comm_mat = readRDS("Data/comm_mat_artiod.RDS"),
                    func_PCoA = readRDS("Data/artiod_func_PCoA.RDS"),
                    iterations = 100)
write.csv(data, "Data/artiod_div_metrics.csv")

# Carnivora
data <- div_metrics(comm_mat = readRDS("Data/comm_mat_carniv.RDS"),
                    func_PCoA = readRDS("Data/carniv_func_PCoA.RDS"),
                    iterations = 100)
write.csv(data, "Data/carniv_div_metrics.csv")

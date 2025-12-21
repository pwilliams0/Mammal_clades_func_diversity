#=================================================================#
# Scatter plots of regional differentiation and diversity metrics #
#=================================================================#

library(tidyverse)
library(sf)
library(usedist)

# Calculate mean phylobetadiversity within and between realms
# Arguments:
#   order = 6-letter taxonomic order abbreviation
#      ("rodent", "chirop", "eulipo", "primat", "carniv", or "artiod")
realm_pb <- function(order) {
 
  # Get cells that contain functional diversity data
  #   These cells have at least as many species as the number of functional axes
  # Also remove cells when there are fewer than 10 cells in the realm
  cell_list <- read.csv(paste("Data/",order,"_div_metrics.csv",sep="")) %>%
    filter(FDis > 0) %>%
    left_join(st_read("Data/cells.shp", quiet=TRUE), "cell_id") %>%
    group_by(Realm) %>%
    mutate(n = n()) %>%
    filter(n > 9) %>%
    # Excludes Neotropical for Eulipotyphla (1 cell)
    # Excludes Panamanian for Primates (4 cells)
    ungroup()
  # Get list of realms where this order occurs
  realm_list <- unique(cell_list$Realm)
  # Load phylobetadiversity distance matrix
  dist <- readRDS(paste("Data/pb_",order,".rds",sep=""))$phylo.beta.sim
  # Keep cells that contain at least one species
  pb_dist <- dist_subset(
    as.dist(dist),
    cell_list$cell_id) %>%
    as.matrix(dimnames=labels(dist))
  # Remove pairwise comparisons within the same cell
  diag(pb_dist) <- NA
  # For each pair of realms, calculate mean phylobetadiversity
  n_realm <- length(realm_list) # number of realms
  realm_pb_mat <- matrix(nrow=n_realm, ncol = n_realm) # matrix to fill in with for loop
  colnames(realm_pb_mat) <- realm_list
  rownames(realm_pb_mat) <- realm_list
  for(i in 1:n_realm){
    for(j in 1:n_realm){
      realm1 <- realm_list[i]
      realm2 <- realm_list[j]
      cell_realm1 <- cell_list %>%
        filter(Realm == realm1) %>%
        pull(cell_id)
      cell_realm2 <- cell_list %>%
        filter(Realm == realm2) %>%
        pull(cell_id)
      # Save mean phylobetadiversity between realms i and j
      realm_pb_mat[i,j] <- mean(na.omit(as.vector(pb_dist[cell_realm1, cell_realm2])))
    }
  }
  # Calculate mean phylobetadiversity within and between realms, and output this
  within_realm <- mean(diag(realm_pb_mat))
  between_realm <- mean(realm_pb_mat[upper.tri(realm_pb_mat, diag=FALSE)])
  output <- setNames(c(within_realm, between_realm),
                     c("within_realm", "between_realm"))
  return(output)
}

# Calculate regional differentiation
#   First calculate mean phylobetadiversity within and between realms for each order
#   Then calculate regional differentiation as:
#      Reg_diff = 1 - (within / between)
Reg_diff <- data.frame(rodent = realm_pb("rodent"),
                   chirop = realm_pb("chirop"),
                   eulipo = realm_pb("eulipo"),
                   primat = realm_pb("primat"),
                   carniv = realm_pb("carniv"),
                   artiod = realm_pb("artiod")) %>%
  rownames_to_column() %>%
  pivot_longer(-rowname) %>%
  pivot_wider(names_from = rowname, values_from = value) %>%
  mutate(Order = name,
         Reg_diff = 1 - (within_realm/between_realm)) %>%
  select(-name)

# Load variance explained by pb only for all orders and each diversity metric
# Join with regional differentiation

# Species richness
SRic <- read.csv("Data/Output/SRic_rodent.csv") %>%
  bind_rows(read.csv("Data/Output/SRic_chirop.csv")) %>%
  bind_rows(read.csv("Data/Output/SRic_eulipo.csv")) %>%
  bind_rows(read.csv("Data/Output/SRic_primat.csv")) %>%
  bind_rows(read.csv("Data/Output/SRic_artiod.csv")) %>%
  bind_rows(read.csv("Data/Output/SRic_carniv.csv")) %>%
  filter(varpart == "Phylobetadiversity only") %>%
  full_join(Reg_diff, by = "Order")

# Funtional richness
FRic <- read.csv("Data/Output/SES.FRic_rodent.csv") %>%
  bind_rows(read.csv("Data/Output/SES.FRic_chirop.csv")) %>%
  bind_rows(read.csv("Data/Output/SES.FRic_eulipo.csv")) %>%
  bind_rows(read.csv("Data/Output/SES.FRic_primat.csv")) %>%
  bind_rows(read.csv("Data/Output/SES.FRic_artiod.csv")) %>%
  bind_rows(read.csv("Data/Output/SES.FRic_carniv.csv")) %>%
  filter(varpart == "Phylobetadiversity only") %>%
  full_join(Reg_diff, by = "Order")

# Functional dispersion
FDis <- read.csv("Data/Output/FDis_rodent.csv") %>%
  bind_rows(read.csv("Data/Output/FDis_chirop.csv")) %>%
  bind_rows(read.csv("Data/Output/FDis_eulipo.csv")) %>%
  bind_rows(read.csv("Data/Output/FDis_primat.csv")) %>%
  bind_rows(read.csv("Data/Output/FDis_artiod.csv")) %>%
  bind_rows(read.csv("Data/Output/FDis_carniv.csv")) %>%
  filter(varpart == "Phylobetadiversity only") %>%
  full_join(Reg_diff, by = "Order")

# Run linear models to check for trends
SR_mod <- lm(R2 ~ Reg_diff, data = SRic)
FR_mod <- lm(R2 ~ Reg_diff, data = FRic)
FDis_mod <- lm(R2 ~ Reg_diff, data = FDis)
summary(SR_mod)
summary(FR_mod)
summary(FDis_mod)

# Plot Reg_diff vs variance explained by pb only for each metric
# Species richness
ggplot(SRic, aes(x=Reg_diff, y=R2, label=Order, color=Order)) +
  theme_classic() + 
  geom_smooth(method="lm", se=FALSE, color="grey50", lwd=.5) +
  geom_point(color=c("black","black","black","red","red","black")) +
  #geom_label() +
  coord_cartesian(ylim=c(.05,.42),
                  xlim=c(.508,.831)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(siz = 12)) +
  labs(x = "Regional differentiation",
       y = "Variance explained by\nphylobetadiversity only")
ggsave("Results/scatter_SRic.svg", width = 3, height = 2.6)
# Functional richness
ggplot(FRic, aes(x=Reg_diff, y=R2, label=Order)) +
  theme_classic() + 
  geom_smooth(method="lm", se=FALSE, color="grey50", lwd=.5) +
  geom_point(color=c("black","red","black","red","red","red")) +
  #geom_label() +
  coord_cartesian(ylim=c(.05,.42),
                  xlim=c(.508,.831)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(siz = 12)) +
  labs(x = "Regional differentiation",
       y = "Variance explained by\nphylobetadiversity only")
ggsave("Results/scatter_FRic.svg", width = 3, height = 2.6)
# Functional dispersion
ggplot(FDis, aes(x=Reg_diff, y=R2, label=Order)) +
  theme_classic() + 
  geom_smooth(method="lm", se=FALSE, color="grey50", lwd=.5) +
  geom_point(color=c("black","black","red","red","red","black")) +
  #geom_label() +
  coord_cartesian(ylim=c(.05,.42),
                  xlim=c(.508,.831)) +
  theme(axis.text = element_text(size = 10),
        axis.title = element_text(siz = 12)) +
  labs(x = "Regional differentiation",
       y = "Variance explained by\nphylobetadiversity only")
ggsave("Results/scatter_FDis.svg", width = 3, height = 2.6)

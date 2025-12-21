#============================#
# Phylobetadiversity figures #
#============================#
library(tidyverse)
library(sf)
library(colorspace)
library(rnaturalearth)
library(rnaturalearthdata)

# ---------- SCREE PLOTS ----------

# Load phylo beta data
# Rodentia
pb_rodent <- readRDS("Data/pb_rodent.rds")$phylo.beta.sim
# Chiroptera
pb_chirop <- readRDS("Data/pb_chirop.rds")$phylo.beta.sim
# Eulipotyphla
pb_eulipo <- readRDS("Data/pb_eulipo.rds")$phylo.beta.sim
# Primates
pb_primat <- readRDS("Data/pb_primat.rds")$phylo.beta.sim
# Carnivora
pb_carniv <- readRDS("Data/pb_carniv.rds")$phylo.beta.sim
# Artiodactyla
pb_artiod <- readRDS("Data/pb_artiod.rds")$phylo.beta.sim

# Calculate stress and make scree plot
# Arguments:
#   pb_dist = phylobetadiversity distance matrix
#   n_axes = max number of axes
#   trymax = max number of runs per set of axes
pb_scree <- function(pb_dist, n_axes, trymax) {
  # Calculate stress for NMDS
  Stress <- c()
  for(i in 1:n_axes){
    set.seed(1234)
    func_mds <- vegan::metaMDS(pb_dist, k=i, trymax=trymax)
    Stress[i] <- func_mds$stress
  }
  # Make scree plot
  par(mar=c(3,3,1,.5))
  plot(1:n_axes, Stress, axes=FALSE, xlab="", ylab="", ylim=c(0,.46))
  abline(h=.2, lty="longdash", col="red")
  lines(1:n_axes, Stress, lwd=.5)
  axis(1, lwd=.5)
  axis(2, lwd=.5)
  title(ylab="Stress", xlab="Number of Dimensions", line=2)
}

# Save scree plots
# Rodentia
svg("Results/PhyloBeta/pb_scree_rodent.svg",
    width=2, height=1.5, pointsize=7)
pb_scree(pb_rodent, 5, 5)
dev.off()
# Chiroptera
svg("Results/PhyloBeta/pb_scree_chirop.svg",
    width=2, height=1.5, pointsize=7)
pb_scree(pb_chirop, 5, 5)
dev.off()
# Eulipotyphla
svg("Results/PhyloBeta/pb_scree_eulipo.svg",
    width=2, height=1.5, pointsize=7)
pb_scree(pb_eulipo, 5, 5)
dev.off()
# Primates
svg("Results/PhyloBeta/pb_scree_primat.svg",
    width=2, height=1.5, pointsize=7)
pb_scree(pb_primat, 5, 5)
dev.off()
# Artiodactyla
svg("Results/PhyloBeta/pb_scree_artiod.svg",
    width=2, height=1.5, pointsize=7)
pb_scree(pb_artiod, 5, 5)
dev.off()
# Carnivora
svg("Results/PhyloBeta/pb_scree_carniv.svg",
    width=2, height=1.5, pointsize=7)
pb_scree(pb_carniv, 5, 5)
dev.off()

# How many axes until stress < 0.2 ?
# Rodent - 3 axes
# Chirop - 3 axes
# Eulipo - 3 axes
# Primat - 1 axis
# Artiod - 3 axes
# Carniv - 2 axes

# ---------- PLOT NMDS AXES ----------

# Set palette for realms
pal_realms <- c('#66CCEE','#ABABAB','#474747','#CCBB44',
               '#228833',"#ABABAB",'#AA3377','#EE6677',
               '#228833','#4477AA',"#EE6677")
# pal_realms <- c("#6929c4","#1192e8","#ee538b","#9f1853",
#                 "#fa4d56","#570408","#198038","#002d9c",
#                 "#005d5d","#b28600","#009d9a")

# Load grid cells
cells <- st_read("Data/cells.shp")

# Plot NMDS axes
# Arguments:
#   nmds = data frame with NMDS coordinates for each cell
#   cells = file with cell id and realm
#   x = NMDS axis to plot on x axis
#   y = NMDS axis to plot on y axis
plot_nmds <- function(nmds, cells, x, y) {
  # Join NMDS with cells to get realm
  pb_plot <- nmds %>% left_join(cells, "cell_id")
  # Define axis limits
  # Add buffer to smaller range so x & y axes have same range
  min_1 <- min(pb_plot[{x}])
  min_2 <- min(pb_plot[{y}])
  max_1 <- max(pb_plot[{x}])
  max_2 <- max(pb_plot[{y}])
  range_1 <- max_1 - min_1
  range_2 <- max_2 - min_2
  range_diff <- abs(range_1 - range_2)
  xlim <- ifelse(range_1 > range_2,
                 list(c(min_1, max_1)),
                 list(c(min_1 - range_diff/2,
                        max_1 + range_diff/2)))[[1]]
  ylim <- ifelse(range_2 > range_1,
                 list(c(min_2, max_2)),
                 list(c(min_2 - range_diff/2,
                        max_2 + range_diff/2)))[[1]]
  # Set palette for realms
  pal <- pal_realms[which(unique(sort(cells$Realm)) %in%
                            unique(sort(pb_plot$Realm)))]
  # Create NMDS plot
  ggplot(pb_plot) +
    aes(x = get({x}), y = get({y}), color = Realm) +
    geom_point(shape = 19, size=.1) +
    scale_color_manual(values=pal) +
    coord_fixed(ratio=1, xlim=xlim, ylim=ylim) +
    theme_classic() +
    labs(x = paste("NMDS", str_sub(x,-1)),
         y = paste("NMDS", str_sub(y,-1))) +
    theme(legend.position = "none",
          axis.text = element_text(color="black",
                                   size=6),
          axis.title = element_text(size=7),
          line = element_line(size = .25))
}

# ----- REALM MAP ----------

# Load coastline for map
world <- ne_coastline(scale = "medium", returnclass = "sf") %>%
  st_transform(crs=st_crs(cells))

# Save map
ggplot() +
  geom_sf(data = world, colour="grey85", linewidth=.1) +
  geom_sf(data = cells, aes(color = Realm, fill = Realm),
          linewidth=.1) +
  scale_color_manual(values=pal_realms) +
  scale_fill_manual(values=pal_realms) +
  theme_void() +
  theme(legend.position = "none")
ggsave("Results/PhyloBeta/realm_map.svg", width = 2.3, height = 1)

# ----- NMDS PLOTS, 2 AXES ----------

# Rodentia
data <- read.csv("Data/pb_NMDS_2_rodent.csv")
plot_nmds(data, cells, "MDS1", "MDS2")
ggsave("Results/PhyloBeta/pb_2axes_rodent.svg", width = 2, height = 2)
# Chiroptera
data <- read.csv("Data/pb_NMDS_2_chirop.csv")
plot_nmds(data, cells, "MDS1", "MDS2")
ggsave("Results/PhyloBeta/pb_2axes_chirop.svg", width = 2, height = 2)
# Eulipotyphla
data <- read.csv("Data/pb_NMDS_2_eulipo.csv") %>%
  mutate(MDS1 = -MDS1)
plot_nmds(data, cells, "MDS1","MDS2")
ggsave("Results/PhyloBeta/pb_2axes_eulipo.svg", width = 2, height = 2)
# Primates
data <- read.csv("Data/pb_NMDS_2_primat.csv") %>%
  mutate(MDS1 = -MDS1)
plot_nmds(data, cells, "MDS1","MDS2")
ggsave("Results/PhyloBeta/pb_2axes_primat.svg", width = 2, height = 2)
# Artiodactyla
data <- read.csv("Data/pb_NMDS_2_artiod.csv")
plot_nmds(data, cells, "MDS1","MDS2")
ggsave("Results/PhyloBeta/pb_2axes_artiod.svg", width = 2, height = 2)
# Carnivora
data <- read.csv("Data/pb_NMDS_2_carniv.csv")
plot_nmds(data, cells, "MDS1","MDS2")
ggsave("Results/PhyloBeta/pb_2axes_carniv.svg", width = 2, height = 2)

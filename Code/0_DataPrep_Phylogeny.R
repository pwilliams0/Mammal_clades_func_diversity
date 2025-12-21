#===============================================#
# Create consensus phylogeny to use in analyses #
#===============================================#

library(ape)
library(tidyverse)
library(phytools)

# Load all Phylacine phylogenies
phylo_all <- ape::read.nexus("Data/Raw/PHYLACINE_1.2/Complete_phylogeny.nex")

# Create consensus tree from all 1000 phylogenies
# Calculate consensus trees of 10 chunks (less computing power)
start <- Sys.time()
mamm_phylo_01 <- consensus.edges(phylo_all[1:100],
                              method = "least.squares")
mamm_phylo_02 <- consensus.edges(phylo_all[101:200],
                              method = "least.squares")
mamm_phylo_03 <- consensus.edges(phylo_all[201:300],
                              method = "least.squares")
mamm_phylo_04 <- consensus.edges(phylo_all[301:400],
                              method = "least.squares")
mamm_phylo_05 <- consensus.edges(phylo_all[401:500],
                              method = "least.squares")
mamm_phylo_06 <- consensus.edges(phylo_all[501:600],
                                 method = "least.squares")
mamm_phylo_07 <- consensus.edges(phylo_all[601:700],
                                 method = "least.squares")
mamm_phylo_08 <- consensus.edges(phylo_all[701:800],
                                 method = "least.squares")
mamm_phylo_09 <- consensus.edges(phylo_all[801:900],
                                 method = "least.squares")
mamm_phylo_10 <- consensus.edges(phylo_all[901:1000],
                                 method = "least.squares")
# Then create consensus tree of the 10 consensus trees
mamm_phylo <- consensus.edges(
  as.multiPhylo(c(mamm_phylo_01, mamm_phylo_02,
                  mamm_phylo_03, mamm_phylo_04,
                  mamm_phylo_05, mamm_phylo_06,
                  mamm_phylo_07, mamm_phylo_08,
                  mamm_phylo_09, mamm_phylo_10)),
  method = "least.squares")
print(Sys.time() - start) # 1.1 hours

# Check that consensus tree is ultrametric
is.ultrametric(mamm_phylo)

# Save consensus phylogeny
write.tree(mamm_phylo, "Data/mamm_phylo.tre")

# Biogeographic isolation leads to functionally divergent communities by restricting the spread of unique traits

### Peter J. Williams, Jedediah F. Brodie, Chia Hsieh, Elise F. Zipkin

---------------------------------

## Abstract

[Final abstract will be added here]

## [Code](Code)

1. **[0_DataPrep_Abiotic.R](Code/0_DataPrep_Abiotic.R)**: Create global grid cells, and calculate environmental variables for each grid cell.
2. **[0_DataPrep_Phylogeny.R](Code/0_DataPrep_Phylogeny.R)**: Create a consensus phylogeny.
3. **[1_DataPrep_Assemblages.R](Code/1_DataPrep_Assemblages.R)**: Create grid cell assemblages for each mammal clade, and process trait data. Some trait data was manually entered after this step (see Tables 2‒6).
4. **[2_DataPrep_DiversityMetrics.R](Code/2_DataPrep_DiversityMetrics.R)**: Calculate species richness, functional richness, and functional dispersion.
5. **[2_DataPrep_PhyloBeta.R](Code/2_DataPrep_PhyloBeta.R)**: Calculate phylobetadiversity.
6. **[3_DataPrep_PhyloBeta_NMDS.R](Code/3_DataPrep_PhyloBeta_NMDS.R)**: Calculate NMDS for phylobetadiversity.
7. **[4_Analysis.R](Code/4_Analysis.R)**: Run variance partitioning analysis.
8. **[5_Figures_BarPlots.R](Code/5_Figures_BarPlots.R)**: Create bar plots of variance partitioning results for Figs. 3a, 4a, & 5a.
9. **[5_Figures_DiversityMaps.R](Code/5_Figures_DiversityMaps.R)**: Create maps of diversity metrics for Figs. 3b, 4b, & 5b. Create maps of mean and SD trait values for Figs. S2‒S7.
10. **[5_Figures_PhyloBeta.R](Code/5_Figures_PhyloBeta.R)**: Create NMDS plots for Fig. 1. Create scree plot for  Fig. S10.
11. **[5_Figures_scatterplots.R](Code/5_Figures_scatterplots.R)**: Calculate regional differentiation, and create scatterplots for Fig. 2.

## [Data](Data)

### Raw data

1. **PHYLACINE_1.2**: PHYLACINE 1.2. Range maps, phylogenies, and trait data for mammals. Avaiable from GitHub ([https://github.com/MegaPast2Future/PHYLACINE_1.2](https://github.com/MegaPast2Future/PHYLACINE_1.2)).
2. **MamFuncDat.txt**: EltonTraits 1.0 data for mammals. Trait data for mammals. Available from Figshare ([https://doi.org/10.6084/m9.figshare.c.3306933.v1](https://doi.org/10.6084/m9.figshare.c.3306933.v1)).
3. **HerbiTraits_1.2.csv**: HerbiTraits 1.2. Trait data for large-bodied terrestrial herbivorous mammals and birds. Available from GitHub ([https://github.com/MegaPast2Future/HerbiTraits](https://github.com/MegaPast2Future/HerbiTraits)).
4. **CarniTraitsv1.5.3.csv**: CarniTraits 1.5.3. Trait data for terrestrial carnivorous mammals. Available from GitHub ([https://github.com/Eamonn-wooster/CarniTraits](https://github.com/Eamonn-wooster/CarniTraits)).
5. **wc2.1_10m_bio**: Bioclimatic variables. Available from WorldClim ([https://www.worldclim.org/data/worldclim21.html](https://www.worldclim.org/data/worldclim21.html)).
6. **mn30_grd**: Global Multi-resolution Terrain Elevation Data (GMTED2010). Available from the United States Geological Survey ([https://topotools.cr.usgs.gov/gmted_viewer/gmted2010_global_grids.php](https://topotools.cr.usgs.gov/gmted_viewer/gmted2010_global_grids.php)).
7. **newRealms**: Zoogeographic realms. Available from the Center for Macroecology, Evolution and Climate at the University of Copenhagen ([https://macroecology.ku.dk/resources/wallace](https://macroecology.ku.dk/resources/wallace)).
8. **ne_50m_land** Land polygons. Available from Natural Earth ([https://www.naturalearthdata.com/downloads/50m-physical-vectors/](https://www.naturalearthdata.com/downloads/50m-physical-vectors/)).

### Derived data products

1. **[cells.shp](Data/cells.shp)**: Grid cells used to define assemblages.
2. **[Landmass_area_cells.csv](Data/Landmass_area_cells.csv)**: Landmass area for each cell. Derived from Williams et al. 2024 landmass area file, available from GitHub ([https://github.com/pwilliams0/Biogeography_and_global_diversity](https://github.com/pwilliams0/Biogeography_and_global_diversity)).
3. **[elev_cells.csv](Data/elev_cells.csv)**: Mean elevation and elevation range of each grid cell.
4. **[climate_cells.csv](Data/climate_cells.csv)**: Climate PCA coordinates for 4 axes for each grid cell.
5. **[mamm_phylo.tre](Data/mamm_phylo.tre)**: Consensus phylogeny.
6. Community matrices for assemblages.
    - **[comm_mat_rodent.RDS](Data/comm_mat_rodent.RDS)**: Rodent community matrix.
    - **[comm_mat_chirop.RDS](Data/comm_mat_chirop.RDS)**: Chiropteran community matrix.
    - **[comm_mat_eulipo.RDS](Data/comm_mat_eulipo.RDS)**: Eulipotyphylan community matrix.
    - **[comm_mat_primat.RDS](Data/comm_mat_primat.RDS)**: Primate community matrix.
    - **[comm_mat_artiod.RDS](Data/comm_mat_artiod.RDS)**: Artiodactyl community matrix.
    - **[comm_mat_carniv.RDS](Data/comm_mat_carniv.RDS)**: Carnivoran community matrix.
7. Trait data for species.
    - **[traits_rodent.csv](Data/traits_rodent.csv)**: Rodent trait data.
    - **[traits_chirop.csv](Data/traits_chirop.csv)**: Chiropteran trait data.
    - **[traits_eulipo.csv](Data/traits_eulipo.csv)**: Eulipotyphylan trait data.
    - **[traits_primat.csv](Data/traits_primat.csv)**: Primate trait data.
    - **[traits_artiod.csv](Data/traits_artiod.csv)**: Artiodactyl trait data.
    - **[traits_carniv.csv](Data/traits_carniv.csv)**: Carnivoran trait data.
    - **[Traits_metadata.pdf](Data/Traits_metadata.pdf)**: Metadata for mammal clades trait data.
8. Coordinates of species in functional space.
    - **[rodent_func_PCoA.RDS](Data/rodent_func_PCoA.RDS)**: Rodent functional space.
    - **[chirop_func_PCoA.RDS](Data/chirop_func_PCoA.RDS)**: Chiropteran functional space.
    - **[eulipo_func_PCoA.RDS](Data/eulipo_func_PCoA.RDS)**: Eulipotyphylan functional space.
    - **[primat_func_PCoA.RDS](Data/primat_func_PCoA.RDS)**: Primate functional space.
    - **[artiod_func_PCoA.RDS](Data/artiod_func_PCoA.RDS)**: Artiodactyl functional space.
    - **[carniv_func_PCoA.RDS](Data/carniv_func_PCoA.RDS)**: Carnivoran functional space.
9. Assemblage-level diversity metrics, inlcuding species richness, functional richness, and functional dispersion.
    - **[rodent_div_metrics.csv](Data/rodent_div_metrics.csv)**: Rodent diversity metrics.
    - **[chirop_div_metrics.csv](Data/chirop_div_metrics.csv)**: Chiropteran diversity metrics.
    - **[eulipo_div_metrics.csv](Data/eulipo_div_metrics.csv)**: Eulipotyphylan diversity metrics.
    - **[primat_div_metrics.csv](Data/primat_div_metrics.csv)**: Primate diversity metrics.
    - **[artiod_div_metrics.csv](Data/artiod_div_metrics.csv)**: Artiodactyl diversity metrics.
    - **[carniv_div_metrics.csv](Data/carniv_div_metrics.csv)**: Carnivoran diversity metrics.
10. Phylobetadiversity data.
    - **[pb_rodent.rds](Data/pb_rodent.rds)**: Rodent phylobetadiversity.
    - **[pb_chirop.rds](Data/pb_chirop.rds)**: Chiropteran phylobetadiversity.
    - **[pb_eulipo.rds](Data/pb_eulipo.rds)**: Eulipotyphylan phylobetadiversity.
    - **[pb_primat.rds](Data/pb_primat.rds)**: Primate phylobetadiversity.
    - **[pb_artiod.rds](Data/pb_artiod.rds)**: Artiodactyl phylobetadiversity.
    - **[pb_carniv.rds](Data/pb_carniv.rds)**: Carnivoran phylobetadiversity.
11. Phylobetadiversity NMDS coordinates using two axes.
    - **[pb_NMDS_2_rodent.csv](Data/pb_NMDS_2_rodent.csv)**: Rodent NMDS coordinates, two axes.
    - **[pb_NMDS_2_chirop.csv](Data/pb_NMDS_2_chirop.csv)**: Chiropteran NMDS coordinates, two axes.
    - **[pb_NMDS_2_eulipo.csv](Data/pb_NMDS_2_eulipo.csv)**: Eulipotyphylan NMDS coordinates, two axes.
    - **[pb_NMDS_2_primat.csv](Data/pb_NMDS_2_primat.csv)**: Primate NMDS coordinates, two axes.
    - **[pb_NMDS_2_artiod.csv](Data/pb_NMDS_2_artiod.csv)**: Artiodactyl NMDS coordinates, two axes.
    - **[pb_NMDS_2_carniv.csv](Data/pb_NMDS_2_carniv.csv)**: Carnivoran NMDS coordinates, two axes.
12. Phylobetadiversity NMDS coordinates using three axes.
    - **[pb_NMDS_3_rodent.csv](Data/pb_NMDS_3_rodent.csv)**: Rodent NMDS coordinates, three axes.
    - **[pb_NMDS_3_chirop.csv](Data/pb_NMDS_3_chirop.csv)**: Chiropteran NMDS coordinates, three axes.
    - **[pb_NMDS_3_eulipo.csv](Data/pb_NMDS_3_eulipo.csv)**: Eulipotyphylan NMDS coordinates, three axes.
    - **[pb_NMDS_3_primat.csv](Data/pb_NMDS_3_primat.csv)**: Primate NMDS coordinates, three axes.
    - **[pb_NMDS_3_artiod.csv](Data/pb_NMDS_3_artiod.csv)**: Artiodactyl NMDS coordinates, three axes.
    - **[pb_NMDS_3_carniv.csv](Data/pb_NMDS_3_carniv.csv)**: Carnivoran NMDS coordinates, three axes.

### [Output](Data/Output)
1. Species richness results from variance partitioning analysis.
    - **[SRic_rodent.csv](Data/Output/SRic_rodent.csv)**: Rodent species richness results.
    - **[SRic_chirop.csv](Data/Output/SRic_chirop.csv)**: Chiropteran species richness results.
    - **[SRic_eulipo.csv](Data/Output/SRic_eulipo.csv)**: Eulipotyphylan species richness results.
    - **[SRic_primat.csv](Data/Output/SRic_primat.csv)**: Primate species richness results.
    - **[SRic_artiod.csv](Data/Output/SRic_artiod.csv)**: Artiodactyl species richness results.
    - **[SRic_carniv.csv](Data/Output/SRic_carniv.csv)**: Carnivoran species richness results.
2. Functional richness results from variance partitioning analysis.
    - **[SES.FRic_rodent.csv](Data/Output/SES.FRic_rodent.csv)**: Rodent functional richness results.
    - **[SES.FRic_chirop.csv](Data/Output/SES.FRic_chirop.csv)**: Chiropteran functional richness results.
    - **[SES.FRic_eulipo.csv](Data/Output/SES.FRic_eulipo.csv)**: Eulipotyphylan functional richness results.
    - **[SES.FRic_primat.csv](Data/Output/SES.FRic_primat.csv)**: Primate functional richness results.
    - **[SES.FRic_artiod.csv](Data/Output/SES.FRic_artiod.csv)**: Artiodactyl functional richness results.
    - **[SES.FRic_carniv.csv](Data/Output/SES.FRic_carniv.csv)**: Carnivoran functional richness results.
2. Functional dispersion results from variance partitioning analysis.
    - **[FDis_rodent.csv](Data/Output/FDis_rodent.csv)**: Rodent functional dispersion results.
    - **[FDis_chirop.csv](Data/Output/FDis_chirop.csv)**: Chiropteran functional dispersion results.
    - **[FDis_eulipo.csv](Data/Output/FDis_eulipo.csv)**: Eulipotyphylan functional dispersion results.
    - **[FDis_primat.csv](Data/Output/FDis_primat.csv)**: Primate functional dispersion results.
    - **[FDis_artiod.csv](Data/Output/FDis_artiod.csv)**: Artiodactyl functional dispersion results.
    - **[FDis_carniv.csv](Data/Output/FDis_carniv.csv)**: Carnivoran functional dispersion results.

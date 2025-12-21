#============================================#
# Bar plots of variance partitioning results #
#============================================#

library(tidyverse)

# Compile results for each order -----------

# Rodentia
results_rodent <- read.csv("Data/Output/SRic_rodent.csv") %>%
  bind_rows(read.csv("Data/Output/SES.FRic_rodent.csv")) %>%
  bind_rows(read.csv("Data/Output/FDis_rodent.csv")) %>%
  mutate(Metric = ifelse (Metric == "SES.FRic", "SES FRic", Metric),
         Metric = factor(Metric, levels = c("SRic","SES FRic","FDis")))
total_R2_rodent <- results_rodent %>%
  group_by(Metric) %>% mutate(total_R2 = sum(R2)) %>% slice(1) %>%
  select(Order, Metric, total_R2)

# Chiroptera
results_chirop <- read.csv("Data/Output/SRic_chirop.csv") %>%
  bind_rows(read.csv("Data/Output/SES.FRic_chirop.csv")) %>%
  bind_rows(read.csv("Data/Output/FDis_chirop.csv")) %>%
  mutate(Metric = ifelse (Metric == "SES.FRic", "SES FRic", Metric),
         Metric = factor(Metric, levels = c("SRic","SES FRic","FDis")))
total_R2_chirop <- results_chirop %>%
  group_by(Metric) %>% mutate(total_R2 = sum(R2)) %>% slice(1) %>%
  select(Order, Metric, total_R2)

# Eulipotyphla
results_eulipo <- read.csv("Data/Output/SRic_eulipo.csv") %>%
  bind_rows(read.csv("Data/Output/SES.FRic_eulipo.csv")) %>%
  bind_rows(read.csv("Data/Output/FDis_eulipo.csv")) %>%
  mutate(Metric = ifelse (Metric == "SES.FRic", "SES FRic", Metric),
         Metric = factor(Metric, levels = c("SRic","SES FRic","FDis")))
total_R2_eulipo <- results_eulipo %>%
  group_by(Metric) %>% mutate(total_R2 = sum(R2)) %>% slice(1) %>%
  select(Order, Metric, total_R2)

# Primates
results_primat <- read.csv("Data/Output/SRic_primat.csv") %>%
  bind_rows(read.csv("Data/Output/SES.FRic_primat.csv")) %>%
  bind_rows(read.csv("Data/Output/FDis_primat.csv")) %>%
  mutate(Metric = ifelse (Metric == "SES.FRic", "SES FRic", Metric),
         Metric = factor(Metric, levels = c("SRic","SES FRic","FDis")))
total_R2_primat <- results_primat %>%
  group_by(Metric) %>% mutate(total_R2 = sum(R2)) %>% slice(1) %>%
  select(Order, Metric, total_R2)

# Artiodactyla
results_artiod <- read.csv("Data/Output/SRic_artiod.csv") %>%
  bind_rows(read.csv("Data/Output/SES.FRic_artiod.csv")) %>%
  bind_rows(read.csv("Data/Output/FDis_artiod.csv")) %>%
  mutate(Metric = ifelse (Metric == "SES.FRic", "SES FRic", Metric),
         Metric = factor(Metric, levels = c("SRic","SES FRic","FDis")))
total_R2_artiod <- results_artiod %>%
  group_by(Metric) %>% mutate(total_R2 = sum(R2)) %>% slice(1) %>%
  select(Order, Metric, total_R2)

# Carnivora
results_carniv <- read.csv("Data/Output/SRic_carniv.csv") %>%
  bind_rows(read.csv("Data/Output/SES.FRic_carniv.csv")) %>%
  bind_rows(read.csv("Data/Output/FDis_carniv.csv")) %>%
  mutate(Metric = ifelse (Metric == "SES.FRic", "SES FRic", Metric),
         Metric = factor(Metric, levels = c("SRic","SES FRic","FDis")))
total_R2_carniv <- results_carniv %>%
  group_by(Metric) %>% mutate(total_R2 = sum(R2)) %>% slice(1) %>%
  select(Order, Metric, total_R2)

# Create barplots -----------

barplot <- function(results){
  data <- results %>% filter(varpart != "Shared") %>%
  mutate(varpart = factor(varpart, levels = c("Phylobetadiversity only",
                                              "Environment only")))
  ggplot(data) +
    aes(x = Metric, y = R2, fill = varpart) +
    geom_col(position=position_dodge(.75), colour="white",
             width=.75, linewidth=0.2) +
    theme_classic() +
    scale_fill_manual(values=c("#8E063B","#023FA5")) +
    scale_y_continuous(limits = c(0,.47), expand = c(0,0)) +
    labs(y = "Variance explained") +
    theme(axis.title.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.text.x=element_text(colour = "black", size=7),
          axis.text.y=element_text(colour = "black", size=7),
          line = element_line(colour = "black", linewidth = 0.2),
          text = element_text(colour = "black", size=7),
          title = element_text(colour = "black", size=7),
          legend.position = "none")
}

# Rodentia
barplot(results_rodent)
ggsave("Results/barplot_rodent.svg", width = 2, height = 2)

# Chiroptera
barplot(results_chirop)
ggsave("Results/barplot_chirop.svg", width = 2, height = 2)

# Eulipotyphla
barplot(results_eulipo)
ggsave("Results/barplot_eulipo.svg", width = 2, height = 2)

# Primates
barplot(results_primat)
ggsave("Results/barplot_primat.svg", width = 2, height = 2)

# Artiodactyla
barplot(results_artiod)
ggsave("Results/barplot_artiod.svg", width = 2, height = 2)

# Carnivora
barplot(results_carniv)
ggsave("Results/barplot_carniv.svg", width = 2, height = 2)

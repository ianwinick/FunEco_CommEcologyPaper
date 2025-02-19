##Ceci Martinez
##February 19 2025

##This is a script to run a preliminary hierarchical cluster analysis to determine 
#functional groups for species occurring in plots within the 
#burn scar of the Museum Fire near Flagstaff, AZ

# load packages -----------------------------------------------------------

library(vegan)
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(dendextend)

# read in data ------------------------------------------------------------
trait_dat <- read_csv("data/TraitTable.csv")
str(trait_dat)

# wrangle data ------------------------------------------------------------

# clean it up, remove last row with NA
trait_dat <- trait_dat %>% 
  slice(-n())

# normalize continuous data
trait_dat_scale <- trait_dat %>%
  mutate(across(c(ldmc, height, sla), ~ as.numeric(scale(.))))

trait_dat_scale <- trait_dat_scale %>% column_to_rownames(var = "spp")

# data exploration --------------------------------------------------------

# some initial plotting/exploration of vars for possible indication of clustering
ggplot(trait_dat, aes(x = sla, y = height, label = spp)) +
  geom_point() +
  geom_text(hjust = 0, nudge_x = 0.09, size = 3) +  
  labs(x = "sla", y = "height") + 
  theme_classic()

ggplot(trait_dat, aes(x = ldmc, y = height, label = spp)) +
  geom_point() +
  geom_text(hjust = 0, nudge_x = 0.003, size = 3) +  
  labs(x = "ldmc", y = "height") + 
  theme_classic()

ggplot(trait_dat, aes(x = ldmc, y = sla, label = spp)) +
  geom_point() +
  geom_text(hjust = 0, nudge_x = 0.003, size = 3) +  
  labs(x = "ldmc", y = "sla") + 
  theme_classic()


###########################################################################
# cluster analysis --------------------------------------------------------
###########################################################################

# create distance matrix based on euclidean distance (can possible explore other options later)

dist_matrix <- dist(trait_dat_scale, method = "euclidean")

# employ different hclust methods (here we use ward, complete linkage, and average linkage methods)
fit_0 <- hclust(dist_matrix, method="ward.D2") #ward method
fit_1 <- hclust(dist_matrix) #complete method 
fit_2 <- hclust(dist_matrix, method = "average") #average method

member_02 <- cutree(fit_0,2)
member_03 <- cutree(fit_0,3)
member_1 <- cutree(fit_1,3)
member_2 <- cutree(fit_2,3)

plot(fit_0) 
plot(fit_1) 
plot(fit_2) 

# visualizing cluster analysis output -----------------------------------------

# plotting dendrogram with 3 gruops
fit0_col3 <- color_branches(as.dendrogram(fit_0), k = 3)
plot(fit0_col)
rect.hclust(fit, k=3, border = "blue")

# plotting dendrogram with 2 gruops
fit0_col2 <- color_branches(as.dendrogram(fit_0), k = 2)
plot(fit0_col)
rect.hclust(fit, k=2, border = "blue")

# plotting clusters on scatterplot
trait_dat$color <- as.factor(member_03)
trait_dat$color2 <- as.factor(member_02)
plot_clusters_sla3 <- ggplot(trait_dat, aes(x = sla, y = ldmc, label = spp, col = color)) + 
  geom_point() + 
  geom_text(vjust = -1, hjust = 0.5, size = 3) +  
  ggtitle("3 Clusters") + 
  theme_classic() + 
  theme(legend.position = "none")
# save ggplot
ggsave("outputs/cluster_scatterplot.pdf", plot_clusters_sla3)

plot_clusters_sla2 <- ggplot(trait_dat, aes(x = sla, y = ldmc, label = spp, col = color2)) + 
  geom_point() + 
  geom_text(vjust = -1, hjust = 0.5, size = 3) +  
  ggtitle("2 Clusters") + 
  theme_classic() 

plot_clusters_ldmc <- ggplot(trait_dat, aes(x = ldmc, y = height, label = spp, col = color)) + 
  geom_point() + 
  geom_text(vjust = -1, hjust = 0.5, size = 3) +  
  theme_classic()


# scree plot for number of groups
scree_plot <- ggplot(fit_0$height %>%
 as.tibble() %>%
    add_column(groups = length(fit_0$height):1) %>%
       rename(height=value),
        aes(x=groups, y=height)) +
        geom_line(col = "cadetblue") + 
        geom_point() +
        ggtitle("ScreePlot") + 
theme_bw() 
# save ggplot
ggsave("outputs/scree_plot.pdf", scree_plot)


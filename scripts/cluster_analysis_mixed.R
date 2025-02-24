##Ceci Martinez
##February 24 2025

##This is a script to run a preliminary hierarchical cluster analysis on both 
#categorical and continuous functonal traits to determine 
#functional groups for species occurring in plots within the 
#burn scar of the Museum Fire near Flagstaff, AZ

# load packages -----------------------------------------------------------

library(vegan)
library(cluster)
library(ggplot2)
library(readr)
library(tidyr)
library(dplyr)
library(tibble)
library(factoextra)
library(dendextend)
library(Rtsne) # for t-SNE plot


# read in data ------------------------------------------------------------

trait_dat <- read_csv("data/TraitTable.csv")
str(trait_dat)

# wrangle data ------------------------------------------------------------

#making relevant groups categorical
trait_dat <- trait_dat %>% 
  mutate(across(c("longevity", "type", "nfix", "resprouting", 
                  "dispersal", "photosynthesis"), as.factor))

# try log transforming first because gower distance is sensitive to non-normality and 
# outliers for continuous variabeles

trait_dat_log <- trait_dat %>%
  column_to_rownames(var = "spp") %>%
  mutate(across(c(ldmc, height, sla), log))
str(trait_dat_log) 

# # normalize continuous data
# trait_dat_log <- trait_dat_log %>%
#   mutate(across(c(ldmc, height, sla), ~ as.numeric(scale(.))))


# data exploration --------------------------------------------------------



###########################################################################
# cluster analysis --------------------------------------------------------
###########################################################################

# create distance matrix based on gower distance (can possible explore other options later)
# https://www.youtube.com/watch?v=2EGgL58drZU

dist_gower <- daisy(trait_dat_log, metric = "gower")

# employ different hclust methods (here we use ward, complete linkage, and average linkage methods)
clust_0 <- hclust(dist_gower, method="ward.D2") #ward method - does best job of giving appx equal cluster size
clust_1 <- hclust(dist_gower, method="single") #single
clust_2 <- hclust(dist_gower, method="complete") #complete

member_02 <- cutree(clust_0,2)
member_03 <- cutree(clust_0,3)
member_04 <- cutree(clust_0,4)

# visualizing cluster analysis --------------------------------------------

# plotting dendrogram with 3 gruops - ward method
fit0_col3 <- color_branches(as.dendrogram(clust_0), k = 3)
plot(fit0_col3)
rect.hclust(clust_0, k=3, border = "blue")

# plotting dendrogram with 4 gropus - ward method
fit0_col4 <- color_branches(as.dendrogram(clust_0), k = 4)
plot(fit0_col4)
rect.hclust(clust_0, k=4, border = "blue")

# plotting dendrogram with 3 gruops - complete method
fit2_col3 <- color_branches(as.dendrogram(clust_2), k = 3)
plot(fit2_col3)
rect.hclust(clust_2, k=3, border = "blue")

# these two different methods yield the same grouping structure

# plotting in reduced dimensions 
tsne_obj <- Rtsne(dist_gower, is_distance = TRUE, perplexity = 6)

tsne_data <- tsne_obj$Y %>%
  data.frame() %>%
  setNames(c("X", "Y")) %>%
  mutate(cluster = factor(member_03),
         name = rownames(trait_dat_log))

viz_ALLtraits <- ggplot(tsne_data, aes(x = X, y = Y)) +
  geom_point(aes(color = factor(member_03)), size = 2.5) +
  geom_text_repel(aes(label = name), size = 3) +
  theme_minimal() +
  scale_color_brewer(palette = "Set1")

ggsave("outputs/cluster_tsneplot_ALLtraits.pdf", viz_ALLtraits)

# adding cluster labels to df
trait_dat_log <- trait_dat_log %>% 
  mutate(cluster = factor(member_03))

# use rf to see which traits are most predictive of clusters
rf_mod <- randomForest(cluster ~ ., data = trait_dat_log, importance = TRUE)
importance(rf_mod)
print(importance(rf_mod)) 
varImpPlot(rf_mod)

# viz cont traits with boxpltos
sla_plot <- ggplot(trait_dat_log, aes(x = cluster, y = sla, fill = cluster)) +
  geom_boxplot() +
  theme_minimal()

ggsave("outputs/sla_plot.pdf", sla_plot)


# viz cat traits --> 
type_plot <- ggplot(trait_dat_log, aes(x = type, fill = cluster)) +
  geom_bar(position = "dodge") +
  theme_minimal()

ggsave("outputs/type_plot.pdf", type_plot)

longevity_plot <- ggplot(trait_dat_log, aes(x = longevity, fill = cluster)) +
  geom_bar(position = "dodge") +
  theme_minimal()

ggsave("outputs/longevity_plot.pdf", longevity_plot)


# scree plot for number of groups
scree_plot <- ggplot(clust_0$height %>%
                       as.tibble() %>%
                       add_column(groups = length(clust_0$height):1) %>%
                       rename(height=value),
                     aes(x=groups, y=height)) +
  geom_line(col = "cadetblue") + 
  geom_point() +
  ggtitle("ScreePlot") + 
  theme_bw() 


# about -------------------------------------------------------------------
# author: maddie wallace
# date: 2-19-25
# run indicator species analysis across burn severity gradient

# packages ----------------------------------------------------------------
#install.packages("indicspecies")
#install.packages("tidyverse")
library(indicspecies)
library(tidyverse)

# analysis ----------------------------------------------------------------
matrix <- read.csv("data/CommunityMatrix.csv")

# reorganize data for pakige
abund <- matrix[,3:ncol(matrix)]
severity <- matrix$Severity

# analyze
indicspp <- multipatt(abund, severity, 
                      func = "r.g", 
                      restcomb = c(1,2,3),
                      control = how(nperm=9999))
summary(indicspp, alpha = 1)

# results
# high severity: VETH, PSMA
# low severity: MUVI
# unburned: nothing with a p < 0.05, but PIPR (p = 0.0753)

#visualize
long_matrix <- matrix |> 
  pivot_longer(cols = 3:ncol(matrix), 
               names_to = "species", 
               values_to = "presence")

ind_spp <- c("VETH", "PSMA", "MUVI", "PIPR")

ind_data <- long_matrix |> 
  filter(species %in% ind_spp)

ggplot(ind_data, aes(x = Severity, y = presence, fill = Severity)) +
  stat_summary(fun = mean, geom = "bar", position = "dodge") +
  facet_wrap(~ species) +
  labs(x = "severity", y = "mean proportion present")



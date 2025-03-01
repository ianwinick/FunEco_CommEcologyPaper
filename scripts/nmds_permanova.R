# Fern Bromley
# Feb 28, 2025
# perform NMDS and do permanova tests on abundance data for museum fire, flagstaff AZ

library(vegan)
library(tidyverse)

abundance <- read_csv("data/CommunityMatrix.csv")

rel_abun <- abundance %>% 
  select(-Severity, -Plot) %>%
  wisconsin()
# rel_abun$Plot <- abundance$Plot
# rel_abun$Severity <- abundance$Severity

plot_info <- select(abundance, Plot, Severity)
severity <- select(abundance, Severity)

# followed this example: https://jonlefcheck.net/2012/10/24/nmds-tutorial-in-r/
# and this: https://uw.pressbooks.pub/appliedmultivariatestatistics/chapter/permanova/

example <- metaMDS(rel_abun, k = 2)

stressplot(example)

burn <- c(rep("Unburned", 20), rep("Low", 18), rep("High", 19))
colors <- c(rep("green", 20), rep("orange", 18), rep("red", 19))

ordiplot(example,type="n")
for(i in unique(burn)) {
  ordihull(example$points[grep(i, burn),], draw="polygon",
           groups=burn[burn==i],col=colors[grep(i,burn)],label=F) }
# orditorp(example,display="species",col="red",air=0.01)
orditorp(example,display="sites",cex=1.25,air=0.01)
legend(col = c("green", "orange", "red"), legend = c("U", "L", "H"),
       x = "topright", pch = 20)

adonis2(rel_abun ~ Severity, data = severity)

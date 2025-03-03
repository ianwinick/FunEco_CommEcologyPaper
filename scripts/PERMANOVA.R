library(tidyverse)
library(vegan)
library(FD)
library(devtools)
install_github("pmartinezarbizu/pairwiseAdonis/pairwiseAdonis")
library(pairwiseAdonis)

data <- read_csv("CommunityMatrix.csv")
traits <- read_csv("TraitTable.csv")%>%
  mutate(ldmc=log(ldmc)) %>%
  mutate(height=log(height)) %>%
  mutate(sla=log(sla)) %>%
  filter(spp %in% colnames(data)) %>%
  column_to_rownames(var="spp")
cwm <- dbFD(traits, comm)$CWM

comm <- data %>%
  select(ARLU:VETH) %>%
  wisconsin()

################################################################################
# TAXONOMIC PERMANOVA ##########################################################
################################################################################

# NMDS stress
tax_nmds <- metaMDS(comm, distance="bray", k=2, trymax=100)
tax_nmds$stress
stressplot(tax_nmds)

# PERMANOVA and post-hoc pairwise PERMANOVA
adonis2(vegdist(comm, method="bray") ~ data$Severity, permutations=9999)
pairwise.adonis(vegdist(comm, method="bray"), data$Severity, perm=9999)

# Test of beta dispersion and post-hoc pairwise test of beta dispersion
anova(betadisper(vegdist(comm, method="bray"), data$Severity, type="centroid"))
TukeyHSD(betadisper(vegdist(comm, method="bray"), data$Severity, type="centroid"))

# Let's get to plotting
tax_scores <- as_tibble(scores(tax_nmds, display="sites"), rownames="sites") %>%
  mutate(severity=case_when(
    sites %in% c(1:20) ~ "Unburned",
    sites %in% c(21:40) ~ "Low",
    sites %in% c(41:60) ~ "High"
  ))
  
tax_centroids <- tax_scores %>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  group_by(severity) %>%
  dplyr::summarize(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))

tax_scores %>%
  mutate(severity=factor(severity, levels=c("Unburned", "Low", "High"))) %>%
  ggplot(aes(x=NMDS1, y=NMDS2, linetype=severity, color=severity)) +
  stat_ellipse(geom="polygon", aes(group=severity, fill=severity), alpha = 0.3, size=.75) +
  geom_point(shape=1, size=2) +
  geom_point(data=tax_centroids, aes(NMDS1, NMDS2), size=3, shape=19) +
  theme_classic()

################################################################################
# FUNCTIONAL PERMANOVA ##########################################################
################################################################################

# Visualize trait distributions 
hist(traits$ldmc)
hist(traits$height)
hist(traits$sla)

# Test if SLA and LDMC are correlated
cor.test(traits$sla, traits$ldmc, method="spearman")
ggplot(traits, aes(sla, ldmc)) +
  geom_point() +
  geom_smooth(method="lm", se=F) +
  theme_classic()

# NMDS stress
fun_nmds <- metaMDS(cwm, distance="euclidean", k=2, trymax=100)
fun_nmds$stress
stressplot(fun_nmds)

# PERMANOVA and post-hoc pairwise PERMANOVA
adonis2(vegdist(cwm, method="euclidean") ~ data$Severity, permutations=9999)
pairwise.adonis(vegdist(cwm, method="euclidean"), data$Severity, perm=9999)

# Test of beta dispersion and post-hoc pairwise test of beta dispersion
anova(betadisper(vegdist(cwm, method="euclidean"), data$Severity, type="centroid"))
TukeyHSD(betadisper(vegdist(cwm, method="euclidean"), data$Severity, type="centroid"))

# Let's get to plotting
fun_scores <- as_tibble(scores(fun_nmds, display="sites"), rownames="sites") %>%
  mutate(severity=case_when(
    sites %in% c(1:20) ~ "Unburned",
    sites %in% c(21:40) ~ "Low",
    sites %in% c(41:60) ~ "High"
  ))
fun_centroids <- fun_scores %>%
  group_by(severity) %>%
  summarise(NMDS1=mean(NMDS1), NMDS2=mean(NMDS2))
fun_scores %>%
  ggplot(aes(x=NMDS1, y=NMDS2, linetype=severity, color=severity)) +
  stat_ellipse(geom="polygon", aes(group=severity, fill=severity), alpha = 0.3, size=.75) +
  geom_point(shape=1, size=2) +
  geom_point(data=fun_centroids, aes(NMDS1, NMDS2), size=3, shape=19) +
  theme_classic()

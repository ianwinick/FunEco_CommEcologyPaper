library(tidyverse)
library(vegan)
library(FD)
library(wesanderson)

data <- read_csv("data/CommunityMatrix.csv")

# traits <- read_csv("data/TraitTable.csv")%>%
#   mutate(ldmc=log(ldmc)) %>%
#   mutate(height=log(height)) %>%
#   mutate(sla=log(sla)) %>%
#   filter(spp %in% colnames(data)) %>%
#   select(spp, sla, height, resprouting, dispersal, seedmass) %>%
#   column_to_rownames(var="spp")
# comm <- data %>%
#   select(ARLU:VETH) %>%
#   decostand(method = "total")

# cwm <- dbFD(traits, comm)$CWM

                  # first approach:
                  # take comm matrix
                  # long form your matrixn (columns plot, spp, abundance (relative))
                  # case_when: species %in% for functional groups
                  # group_by severity, functional group
                  # summarise to average the 
                  
                  # comm$severity <- data$Severity
                  # group_cover <- comm %>% 
                  #   pivot_longer(cols = ARLU:VETH, names_to = "spp") %>% 
                  #   mutate(group = case_when(spp %in% c("CARO", "ELEL", "FEAR",
                  #                                       "MUMO", "MUVI", "PIPR", 
                  #                                       "SCSC") ~ "graminoid",
                  #                            spp == "CEFE" ~ "shrub",
                  #                            spp == "QUGA" ~ "tree",
                  #                            .default = "forb")) %>% 
                  #   group_by(severity, group) %>% 
                  #   summarise(mean_cov = mean(value))
                  # 
                  # group_cover$severity <- factor(group_cover$severity, c("U", "L", "H"))
                  # 
                  # ggplot(group_cover, aes(x = severity, y = mean_cov*100, fill = group))+
                  #   geom_bar(stat = "identity")+
                  #   labs(x = "Burn severity", y = "Relative Cover (%)",
                  #        fill = "Functional type")+
                  #   theme_light()

# ian's second suggested approach:

# use absolute cover
# assign fxnl groups
# sum functional groups within each plot
# average across plots within severities
# THEN relativize the averaged plot by severity

# my thoughts:
# i think we should sum absolute functional groups within each plot AND 
# then relativize the functional groups within severities so it adds to 100%

comm <- data
# comm$severity <- data$Severity
# comm$plot <- data$Plot

group_cover2 <- comm %>% 
  pivot_longer(cols = ARLU:VETH, names_to = "spp") %>% 
  mutate(group = case_when(spp %in% c("CARO", "ELEL", "FEAR",
                                      "MUMO", "MUVI", "PIPR", 
                                      "SCSC") ~ "Graminoid",
                           spp == "CEFE" ~ "Shrub",
                           spp == "QUGA" ~ "Tree",
                           .default = "Forb")) %>% 
  group_by(Severity, Plot, group) %>% 
  summarise(fun_cov = sum(value)) %>% 
  ungroup()

fun_cover <- group_cover2 %>% 
  pivot_wider(names_from = group, values_from = fun_cov)

rel_fun_cover <- fun_cover %>% 
  select(-Severity, -Plot) %>% 
  decostand(method = "total")

rel_fun_cover$severity <- fun_cover$Severity

cover_df <- rel_fun_cover %>% 
  pivot_longer(cols = Forb:Tree, names_to = "group") %>% 
  group_by(severity, group) %>% 
  summarise(cov = mean(value))
cover_df$severity <- factor(cover_df$severity, c("U", "L", "H"))
# cover_df$severity <- factor(cover_df$severity, c("U", "L", "H"))


# pal <- wes_palette("Rushmore1", 4, type = "discrete")
  
ggplot(cover_df, aes(x = severity, y = cov*100, fill = group))+
  theme_light(base_size = 26)+
  geom_bar(stat = "identity", color = "black")+
  labs(x = "Burn severity", y = "Relative Cover (%)",
       fill = "Functional type")+
  scale_fill_manual(values = wes_palette("Cavalcanti1", 4))

#group_cover$severity <- factor(group_cover$severity, c("U", "L", "H"))








# attempt 3: DON'T USE
comm <- data

group_cover3 <- comm %>% 
  pivot_longer(cols = ARLU:VETH, names_to = "spp") %>% 
  mutate(group = case_when(spp %in% c("CARO", "ELEL", "FEAR",
                                      "MUMO", "MUVI", "PIPR", 
                                      "SCSC") ~ "graminoid",
                           spp == "CEFE" ~ "shrub",
                           spp == "QUGA" ~ "tree",
                           .default = "forb")) %>% 
  group_by(Severity, Plot, group) %>% 
  summarise(fun_cov = sum(value)) %>% 
  ungroup()

avg_fun_cover <- group_cover3 %>% 
  group_by(Severity, group) %>% 
  summarise(mean_cov = mean(fun_cov)) %>% 
  ungroup()

rel_fun_cover <- avg_fun_cover %>% 
  pivot_wider(names_from = group, values_from = mean_cov) %>% 
  select(-Severity) %>% 
  wisconsin()
  
rel_fun_cover$severity <- c("H", "L", "U")

cover_df <- rel_fun_cover %>%
  pivot_longer(cols = forb:tree, names_to = "group")
cover_df$severity <- factor(cover_df$severity, c("U", "L", "H"))

ggplot(cover_df, aes(x = severity, y = value*100, fill = group))+
  geom_bar(stat = "identity")+
  labs(x = "Burn severity", y = "Relative Cover (%)",
       fill = "Functional type")+
  scale_fill_grey()+
  theme_light()




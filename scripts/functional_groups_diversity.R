library(tidyverse)

data <- read_csv("data/CommunityMatrix.csv")

traits <- read_csv("data/TraitTable.csv")
comm <- data %>%
  select(ARLU:VETH) %>%
  wisconsin()
cwm <- dbFD(traits, comm)$CWM

# first approach:
# take comm matrix
# long form your matrixn (columns plot, spp, abundance (relative))
# case_when: species %in% for functional groups
# group_by severity, functional group
# summarise to average the 

comm$severity <- data$Severity
group_cover <- comm %>% 
  pivot_longer(cols = ARLU:VETH, names_to = "spp") %>% 
  mutate(group = case_when(spp %in% c("CARO", "ELEL", "FEAR",
                                      "MUMO", "MUVI", "PIPR", 
                                      "SCSC") ~ "graminoid",
                           spp == "CEFE" ~ "shrub",
                           spp == "QUGA" ~ "tree",
                           .default = "forb")) %>% 
  group_by(severity, group) %>% 
  summarise(mean_cov = mean(value))

group_cover$severity <- factor(group_cover$severity, c("U", "L", "H"))

ggplot(group_cover, aes(x = severity, y = mean_cov*100, fill = group))+
  geom_bar(stat = "identity")+
  labs(x = "Burn severity", y = "Relative Cover (%)",
       fill = "Functional type")+
  theme_light()

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
                                      "SCSC") ~ "graminoid",
                           spp == "CEFE" ~ "shrub",
                           spp == "QUGA" ~ "tree",
                           .default = "forb")) %>% 
  group_by(Severity, Plot, group) %>% 
  summarise(fun_cov = sum(value)) %>% 
  ungroup()

rel_fun_cover <- group_cover2 %>% 
  pivot_wider(names_from = group, values_from = fun_cov) %>% 
  select(-Severity, -Plot) %>% 
  wisconsin()
rel_fun_cover$severity <- data$Severity

cover_df <- rel_fun_cover %>% 
  pivot_longer(cols = forb:tree, names_to = "group") %>% 
  group_by(severity, group) %>% 
  summarise(cov = mean(value))
cover_df$severity <- factor(cover_df$severity, c("U", "L", "H"))

ggplot(cover_df, aes(x = severity, y = cov*100, fill = group))+
  geom_bar(stat = "identity")+
  labs(x = "Burn severity", y = "Relative Cover (%)",
       fill = "Functional type")+
  scale_fill_grey()+
  theme_light()

#group_cover$severity <- factor(group_cover$severity, c("U", "L", "H"))

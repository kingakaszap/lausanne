library(readxl)
library(tidyverse)
diversity <- read_excel("applied_ecology/diversity.xlsx")
View(diversity)

# species richness

diversity_long <- diversity %>% 
  pivot_longer(cols = 2:19, names_to = "Site", values_to = "Presence")
View(diversity_long)
diversity_long <- na.omit(diversity_long)

(richness <-  diversity_long %>% 
  group_by(Site) %>% 
  summarise(sp.richness = length(unique(Species))))

(barplot_richness <- ggplot(richness, aes(x = reorder(Site, sp.richness), y = sp.richness, fill = Site)) +
    # setting the x and y axes, and asking R to use different colors for each park
    geom_bar (stat = "identity") +
    # specifying that we want a barplot ?
    labs (x = "\nPark", y = "Species richness\n") +
    # giving informative names to the axes
    theme_classic() +
    #adding a color palette - feel free to add a different one!
    theme (legend.position = "none",
           # removing the legend as the axes provide enough information 
           axis.text = element_text(size = 14), 
           axis.title = element_text(size = 16),
           #increasing font size
           axis.text.x = element_text(angle = 45, hjust = 1)))
#tilting the text of the x axis

ggsave(barplot_richness, file = "applied_ecology/richnessplot.png", width = 6, height = 6)

# shannons diversity
diversity_long <- diversity_long %>% 
  group_by(Site) %>% 
  mutate(relative.abundance = Presence/sum(Presence)) 

shannon <- diversity_long %>% 
  group_by(Site) %>% 
  summarise(shannons.div = -sum(relative.abundance*log(relative.abundance))) %>%
  ungroup()

(barplot_shannon <- ggplot(shannon, aes(x = reorder(Site, shannons.div), y = shannons.div, fill = Site)) +
    # setting the x and y axes, and asking R to use different colors for each park
    geom_bar (stat = "identity") +
    # specifying that we want a barplot ?
    labs (x = "\nSite", y = "Species richness\n") +
    # giving informative names to the axes
    theme_classic() +
    #adding a color palette - feel free to add a different one!
    theme (legend.position = "none",
           # removing the legend as the axes provide enough information 
           axis.text = element_text(size = 14), 
           axis.title = element_text(size = 16),
           #increasing font size
           axis.text.x = element_text(angle = 45, hjust = 1)))
#tilting the text of the x axis

# Calculate H'max (maximum possible Shannon diversity for each site)
richness <- richness %>% 
  mutate(H_max = ifelse(sp.richness > 1, log(sp.richness), NA))  # Avoid log(1) issue

# Merge Shannon diversity (H') with richness (S) per site
evenness <- left_join(shannon, richness, by = "Site") %>%
  mutate(Pielou_Evenness = shannons.div / H_max) %>%
  select(Site, Pielou_Evenness)  # Keep only relevant columns

# Print results
print(evenness)

# Plot Pielou's Evenness per site
(barplot_evenness <- ggplot(evenness, aes(x = reorder(Site, Pielou_Evenness), y = Pielou_Evenness, fill = Site)) +
    geom_bar(stat = "identity") +
    labs(x = "\nSite", y = "Pielou's Evenness (J')\n") +
    theme_classic() +
    theme(legend.position = "none",
          axis.text = element_text(size = 14), 
          axis.title = element_text(size = 16),
          axis.text.x = element_text(angle = 45, hjust = 1)))

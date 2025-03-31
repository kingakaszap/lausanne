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

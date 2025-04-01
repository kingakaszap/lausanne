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

max(richness$sp.richness)

(barplot_richness <- ggplot(richness, aes(x = reorder(Site, sp.richness), y = sp.richness, fill = Site)) +
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
           axis.title.x = element_blank(),
           #increasing font size
           axis.text.x = element_text(angle = 45, hjust = 1)))
#tilting the text of the x axis

ggsave(barplot_richness, file = "applied_ecology/richnessplot.png", width = 8, height = 6)

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
    labs (x = "\nSite", y = "Shannon's diversity\n") +
    # giving informative names to the axes
    theme_classic() +
    #adding a color palette - feel free to add a different one!
    theme (legend.position = "none",
           # removing the legend as the axes provide enough information 
           axis.text = element_text(size = 14), 
           axis.title = element_text(size = 16),
           axis.title.x = element_blank(),
           #increasing font size
           axis.text.x = element_text(angle = 45, hjust = 1)))
#tilting the text of the x axis

ggsave(barplot_shannon, file = "applied_ecology/shannonplot.png", width = 8, height = 6)


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
          axis.title.x= element_blank(),
          axis.text = element_text(size = 14), 
          axis.title = element_text(size = 16),
        axis.text.x = element_text(angle = 45, hjust = 1)))
ggsave(barplot_evenness, file = "applied_ecology/evennessplot.png", width = 8, height = 6)



# with vegan
library(vegan)
df <- diversity
# Extract site names and remove the first column (species IDs)
site_names <- df[1, -1]  # Assuming the first row has site names
species_data <- df[-1, ]  # Remove first row (site names)

# Convert species data to numeric, replacing NAs with 0
species_data <- species_data %>% mutate_all(~as.numeric(replace(., is.na(.), 0)))

# Transpose data so that sites are rows and species are columns
colnames(species_data) <- site_names
species_data <- as.data.frame(t(species_data))
# Calculate species richness (number of species per site)
richness <- specnumber(species_data)

# Calculate Shannon's diversity index
shannon <- diversity(species_data, index = "shannon")

# Calculate Pielou's evenness (Shannon's index divided by log species richness)
evenness <- ifelse(richness > 1, shannon / log(richness), NA)

# Combine results into a dataframe
results <- data.frame(Site = rownames(species_data),
                      Richness = richness,
                      Shannon = shannon,
                      Evenness = evenness)

# View results
print(results)

rankabundance <- diversity_long %>% 
  # making new dataframe for rank abundance plots
  group_by (Site) %>% 
  mutate(rel.abund.percent = relative.abundance * 100,
         # adding new column for % relative abundance  
         rank = (rank(- Presence, ties.method="random"))) %>% 
  # new columnn for rank  
  ungroup()
# remove groupings

(rank_abundance_plots <- rankabundance %>% 
    ggplot (aes (x = rank, y = rel.abund.percent)) +
    # specifying the x and y axes
    geom_point (aes (color = Species, fill = Species), size = 2) +
    # specifying that we want a scatter plot
    geom_line (colour= "black") +
    # adding a line connecting the points
    labs (x = "\nRank", y = "Relative abundance (%)\n") +
    facet_wrap (~Site) +
    # making separate plots for sites
    theme_classic() +
    # adding a theme from ggthemes - feel free to add a different one!
    theme (legend.position = "none"))
# removing the legend as the axes provide enough information

ggsave(rank_abundance_plots, file = "applied_ecology/rankabundance.png", dpi = 600)

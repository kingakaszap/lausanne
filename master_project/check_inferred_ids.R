# compare inferred id-s against id-s in vcf (there shouldnt be)
# and when they were born

# libraries and data ----
library(tidyverse)
library(readr)
library(lubridate)


pedigree <- read_csv("master_project/data/coreccted_consensus_pedigree.txt")
# add FamilyID column to pedigree based on parents
pedigree$FamilyID <- paste(pedigree$MumId, pedigree$DadId, sep = "_")

sequenced <- read_csv("master_project/data/RingIds_seq1.txt") # individuals who are sequenced
in_vcf <- read.table("master_project/data/ids_in_vcf.txt", header = FALSE, sep = "")
in_vcf<- in_vcf %>% rename(RingId = V1)
inferred <- read_csv("master_project/data/ids_inferred_snipar.csv")
inferred<- as.data.frame(inferred)
str(inferred)
inferred<- inferred %>% rename( RingId =  ID)
length(intersect(inferred$ID, sequenced)) # SUPER

bird <- read.csv("master_project/data/bird.csv")
str(bird$HatchDate)

bird$HatchDate <- dmy_hms(bird$HatchDate)   
bird$HatchDate <- as.Date(bird$HatchDate) # strip time

hatchdate_df <- bird %>% select(BirdId, RingId, HatchDate, LastObservedAlive )

inferred <- left_join(inferred, hatchdate_df, by = "RingId")
write_csv(inferred,"master_project/data/inferred_with_years.csv")


# plot hatch years
inferred <- inferred %>%
  mutate(hatchyear = year(HatchDate))

(which_years <- ggplot(inferred, aes(x = hatchyear)) +
  geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
  labs(x = "year hatched",y = "number of inferred individuals") +
  theme_classic())
ggsave("master_project/plots/which_year_inferred_individuals.png", dpi = 600)

# for id-s in vcf

in_vcf <- left_join(in_vcf, hatchdate_df, by = "RingId")
in_vcf <- in_vcf %>%
  mutate(hatchyear = year(HatchDate))
(which_years_vcf <- ggplot(in_vcf, aes(x = hatchyear)) +
    geom_histogram(binwidth = 1, fill = "steelblue", color = "white") +
    labs(x = "year hatched",y = "number of individuals") +
    theme_classic())
ggsave("master_project/plots/which_year_vcf_individuals.png", dpi = 600)



# chatgpt attempt - make 1 plot

# Prepare inferred dataset
inferred_plot <- inferred %>% mutate(source = "Inferred") 

# Prepare VCF dataset
vcf_plot <- in_vcf %>% mutate(source = "Sequenced")  

combined_df <- bind_rows(inferred_plot, vcf_plot)

# Plot
(combined_hist <- ggplot(combined_df, aes(x = hatchyear, fill = source)) +
    geom_histogram(data = subset(combined_df, source == "Inferred"),
                   binwidth = 1, fill = "steelblue", alpha = 1, color = "white") +
    geom_histogram(data = subset(combined_df, source != "Inferred"),
                   binwidth = 1, fill = "grey70", alpha = 0.5, color = "white") +
    labs(x = "Year hatched", y = "Number of individuals", fill = "Dataset") +
    theme_classic())

ggsave("master_project/plots/which_year_combined.png", combined_hist, dpi = 600)


# find families that have 2 or more genotyped children and both parents genotyped
# after -> find grandparents etc like chains

library(tidyverse)
library(readr)
library(readr)
pedigree <- read_csv("master_project/data/coreccted_consensus_pedigree.txt")
View(pedigree)

sequenced <- read_csv("master_project/data/RingIds_seq1.txt")
View(sequenced)

full_families <- pedigree %>%
  filter(!is.na(MumId) & !is.na(DadId))


# Convert IDs to numeric (fixes any character/numeric mismatch issues)
pedigree <- pedigree %>%
  mutate(across(c(RingId, MumId, DadId), as.numeric))

# Find families where both parents are known
valid_families <- pedigree %>%
  filter(!is.na(MumId) & !is.na(DadId))








# Find sequenced children
sequenced_families <- valid_families %>%
  inner_join(pedigree, by = c("MumId", "DadId")) %>%  # Join to find children
  filter(RingId.y %in% sequenced$RingId) %>%  # Keep only sequenced children
  group_by(MumId, DadId) %>%
  summarise(Children = list(RingId.y), .groups = "drop") %>%  # Collect children in a list
  filter(lengths(Children) >= 2)  # Keep only families with at least 2 sequenced children

# Ensure there are sequenced families before proceeding
if (nrow(sequenced_families) > 0) {
  max_children <- max(lengths(sequenced_families$Children)) 
  
  output <- sequenced_families %>%
    mutate(Children = lapply(Children, function(x) c(x, rep(NA, max_children - length(x))))) %>%
    unnest_wider(Children, names_sep = "_")  # Spread children into separate columns
  
  print(output)
} else {
  print("didn't work")
}
print(pedigree$RingId[pedigree$RingId %in% sequenced])
str
View(output)



family_list <- full_families %>%
  rowwise() %>%
  mutate(Children = list(RingId[RingId %in% sequenced])) %>%
  ungroup()

filtered_families <- family_list %>%
  filter(length(Children) >= 2) %>%
  select(MumId, DadId, Children)


max_children <- max(sapply(filtered_families$Children, length))
families <- filtered_families %>%
  mutate(Children = lapply(Children, function(x) c(x, rep(NA, max_children - length(x))))) %>%
  tidyr::unnest_wider(Children, names_sep = "_")

View(families)
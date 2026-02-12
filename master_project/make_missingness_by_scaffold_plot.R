# combine all the scaffolds and get the by-position missingness

# libraries ----
.libPaths(c("/users/kkaszap/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(tidyverse)

root_dir <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/na_summaries"

# list all csv files that match "na_by_snp_" recursively
snp_files <- list.files(path = root_dir,pattern = "^na_by_snp_.*\\.csv$",recursive = TRUE,
  full.names = TRUE)

# read and combine
all_snps <- snp_files %>%lapply(read.csv) %>% bind_rows()

# make sure position is numeric
all_snps$position <- as.numeric(as.character(all_snps$position))
head(all_snps)

scaffold_order <- unique(all_snps$scaffold_id)
all_snps$scaffold_id <- factor(all_snps$scaffold_id, levels = scaffold_order)

# create simple continuous x-axis for plotting
all_snps <- all_snps %>%
  group_by(scaffold_id) %>%
  arrange(position) %>%
  mutate(pos_cum = row_number()) %>%
  ungroup()

# midpoints for scaffold labels
scaffold_mid <- all_snps %>%
  group_by(scaffold_id) %>%
  summarize(mid = mean(pos_cum))

# -----------------------
# 1. Faceted plot
# -----------------------
facet_plot <- ggplot(all_snps, aes(x = position, y = percent_missing_in_individuals)) +
  geom_point(aes(color = percent_missing_in_individuals == 0), 
             size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue"), guide = "none") +
  facet_wrap(~ scaffold_id, scales = "free_x") +
  theme_classic() +
  labs(x = "SNP position", y = "% of individuals missing",
       title = "Missingness by SNP position (faceted by scaffold)") +
  scale_y_continuous(limits = c(0, 100))

# save faceted plot
ggsave(file.path(root_dir, "missingness_faceted_by_scaffold.png"), facet_plot, dpi = 300, width = 12, height = 6)

# -----------------------
# 2. Combined plot
# -----------------------

# order scaffolds numerically by extracting number
all_snps <- all_snps %>%
  arrange(scaffold_id, position) %>%
  group_by(scaffold_id) %>%
  mutate(pos_cum = row_number()) %>%
  ungroup()

# Offset cumulative positions so scaffolds are continuous
scaffold_max <- all_snps %>%
  group_by(scaffold_id) %>%
  summarise(max_pos = max(pos_cum), .groups = "drop") %>%
  mutate(cum_offset = lag(cumsum(max_pos), default = 0))

all_snps <- all_snps %>%
  left_join(scaffold_max %>% select(scaffold_id, cum_offset), by = "scaffold_id") %>%
  mutate(pos_cum = pos_cum + cum_offset)

# Get midpoints for scaffold labels
scaffold_ticks <- all_snps %>%
  group_by(scaffold_id) %>%
  summarise(mid = mean(pos_cum), .groups = "drop")

# Plot combined
combined_plot <- ggplot(all_snps, aes(x = pos_cum, y = percent_missing_in_individuals)) +
  geom_point(aes(color = percent_missing_in_individuals == 0), 
             size = 0.6, alpha = 0.6) +
  scale_color_manual(values = c("FALSE" = "black", "TRUE" = "blue"), guide = "none") +
  scale_x_continuous(breaks = scaffold_ticks$mid, labels = scaffold_ticks$scaffold_id) +
  theme_classic() +
  labs(x = "Scaffold", y = "% of individuals missing",
       title = "Missingness by SNP position across all scaffolds") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_y_continuous(limits = c(0, 100))

# save combined plot
ggsave(file.path(root_dir, "missingness_combined_all_scaffolds.png"), combined_plot, dpi = 300, width = 14, height = 6)
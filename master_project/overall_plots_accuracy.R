.libPaths(c("/users/kkaszap/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(tidyverse)

# ---------- paths ----------
scaffold_summary_dir <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/accuracy_summaries"
out_plot_dir <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/summary_plots"
if (!dir.exists(out_plot_dir)) dir.create(out_plot_dir, recursive = TRUE)

# ---------- 1. Mean accuracy per scaffold ----------
scaffold_files <- list.files(scaffold_summary_dir, pattern = "scaffold_summary_.*\\.csv$", full.names = TRUE)

scaffold_summary_all <- lapply(scaffold_files, read.csv) %>%
  bind_rows()

# Plot mean absolute difference per scaffold
plot_scaffold <- ggplot(scaffold_summary_all, aes(x = scaffold_id, y = mean_abs_dif)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_classic() +
  labs(title = "Mean Accuracy per Scaffold",
       x = "Scaffold", y = "Mean Absolute Difference") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(file.path(out_plot_dir, "mean_accuracy_per_scaffold.png"),
       plot = plot_scaffold, dpi = 300, width = 10, height = 4)

# ---------- 2. Genome-wide SNP accuracy ----------
snp_files <- list.files(scaffold_summary_dir, pattern = "accuracy_by_snp_.*\\.csv$", full.names = TRUE)

snp_all <- lapply(snp_files, function(f) {
  df <- read.csv(f)
  df$scaffold_id <- sub("accuracy_by_snp_|\\.csv", "", basename(f))
  df
}) %>% bind_rows()

plot_snp_genome <- ggplot(snp_all, aes(x = SNP, y = mean_abs_dif, color = scaffold_id)) +
  geom_point(size = 0.5, alpha = 0.6) +
  theme_classic() +
  labs(title = "Genome-wide SNP Accuracy",
       x = "SNP position", y = "Mean Absolute Difference",
       color = "Scaffold") +
  theme(axis.text.x = element_blank())
ggsave(file.path(out_plot_dir, "genome_wide_snp_accuracy.png"),
       plot = plot_snp_genome, dpi = 300, width = 12, height = 4)

# ---------- 3. Mean accuracy per individual ----------
ind_files <- list.files(scaffold_summary_dir, pattern = "accuracy_by_individual_.*\\.csv$", full.names = TRUE)

ind_all <- lapply(ind_files, read.csv) %>%
  bind_rows() %>%
  group_by(Sample) %>%
  summarise(mean_abs_dif = mean(mean_abs_dif, na.rm = TRUE), .groups = "drop") %>%
  filter(!is.na(mean_abs_dif))  # remove individuals with all NA

plot_ind <- ggplot(ind_all, aes(x = Sample, y = mean_abs_dif)) +
  geom_bar(stat = "identity", fill = "coral") +
  theme_classic() +
  labs(title = "Mean Accuracy per Individual Across Scaffolds",
       x = "Sample", y = "Mean Absolute Difference") +
  theme(axis.text.x = element_blank())
ggsave(file.path(out_plot_dir, "mean_accuracy_per_individual.png"),
       plot = plot_ind, dpi = 300, width = 10, height = 4)

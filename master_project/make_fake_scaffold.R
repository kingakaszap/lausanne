# fake accuracy comparison
# what if we always predict 2

.libPaths(c("/users/kkaszap/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(tidyverse)

accuracy_root <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/accuracy_summaries"

# -------------------------------------------------
# Find all comparison_long_df files
# -------------------------------------------------
comparison_files <- list.files(
  accuracy_root,
  pattern = "^comparison_long_df_.*\\.csv$",
  recursive = TRUE,
  full.names = TRUE)

stopifnot(length(comparison_files) > 0)

# -------------------------------------------------
# Compute fake mean abs error per scaffold
# -------------------------------------------------
fake_error_by_scaffold <- lapply(comparison_files, function(f) {
  
  df <- read.csv(f)
  
  tibble(
    scaffold_id = unique(df$scaffold_id),
    fake_mean_abs_dif = mean(abs(2 - df$Actual), na.rm = TRUE) )}) %>%
  bind_rows() %>%
  arrange(scaffold_id)

# Save
write.csv(
  fake_error_by_scaffold,
  file.path(accuracy_root, "fake_mean_abs_error_all_scaffolds.csv"),
  row.names = FALSE)

fake_error_by_scaffold

real_summary <- read.csv(file.path(
  accuracy_root, "all_scaffolds_summary.csv"))

comparison <- left_join(
  real_summary,
  fake_error_by_scaffold,
  by = "scaffold_id")

comparison_long <- comparison %>%
  select(scaffold_id, mean_abs_dif, fake_mean_abs_dif) %>%
  pivot_longer(
    cols = c(mean_abs_dif, fake_mean_abs_dif),
    names_to = "type",
    values_to = "mean_abs_dif")

plot_fake_vs_real <- ggplot(comparison_long,
       aes(x = scaffold_id, y = mean_abs_dif, fill = type)) +
  geom_bar(stat = "identity", position = "dodge") +
  theme_classic() +
  labs(
    title = "Real vs majority (always = 2) mean error",
    x = "Scaffold",
    y = "Mean absolute difference") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

ggsave(filename = file.path(accuracy_root, "barplot_fake_2_s_scaffold.png"),
       plot = plot_fake_vs_real,
       dpi = 300, width = 12, height = 5)
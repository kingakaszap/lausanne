# -----------------------------------------
# Post-processing: Combine scaffold summaries
# -----------------------------------------
# libraries ----
.libPaths(c("/users/kkaszap/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(rhdf5)
library(SNPRelate) # remember to install it on cluster with an interactive session ! 
library(tidyverse) # same

# Folder where per-scaffold summaries are saved
summary_folder <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/accuracy_summaries"

# Output files
combined_csv <- file.path(summary_folder, "all_scaffolds_summary.csv")
barplot_png  <- file.path(summary_folder, "all_scaffolds_mean_abs_diff.png")

# -----------------------------------------
# Step 1: list all scaffold summary CSVs
# -----------------------------------------
summary_files <- list.files(summary_folder, pattern = "^scaffold_summary_.*\\.csv$", full.names = TRUE)

if(length(summary_files) == 0) stop("No scaffold summary CSVs found in ", summary_folder)

# -----------------------------------------
# Step 2: read and combine
# -----------------------------------------
all_scaffolds_summary <- summary_files %>%
  lapply(read.csv, stringsAsFactors = FALSE) %>%
  bind_rows() %>%
  arrange(scaffold_id)

# Optional: make scaffold_id a factor to preserve order in plots
all_scaffolds_summary$scaffold_id <- factor(all_scaffolds_summary$scaffold_id,
                                            levels = unique(all_scaffolds_summary$scaffold_id))

# -----------------------------------------
# Step 3: save combined CSV
# -----------------------------------------
write.csv(all_scaffolds_summary, file = combined_csv, row.names = FALSE)
cat("Saved combined scaffold summary CSV to:", combined_csv, "\n")

# -----------------------------------------
# Step 4: Bar plot of mean_abs_dif per scaffold
# -----------------------------------------
plot_scaffold_accuracy <- ggplot(all_scaffolds_summary, aes(x = scaffold_id, y = mean_abs_dif)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  theme_classic() +
  labs(title = "Mean absolute imputation error per scaffold",
       x = "Scaffold",
       y = "Mean absolute difference") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# Save plot
ggsave(filename = barplot_png,
       plot = plot_scaffold_accuracy,
       dpi = 300, width = 12, height = 5)

cat("Saved bar plot of mean_abs_dif per scaffold to:", barplot_png, "\n")

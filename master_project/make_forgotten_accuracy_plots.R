# do forgotten accuracy - by - snp plots
# 2026 jan–feb
.libPaths(c("/users/kkaszap/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(rhdf5)
library(tidyverse)

#  paths ----------
accuracy_dir <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/accuracy_summaries"

# get SLURM array id ----------
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) stop("No SLURM_ARRAY_TASK_ID found")

#  find scaffold folders ----------
scaffold_dirs <- list.dirs(accuracy_dir, recursive = FALSE, full.names = TRUE)

if (task_id > length(scaffold_dirs)) {stop("Task ID exceeds number of scaffolds")}

scaffold_dir <- scaffold_dirs[task_id]
scaffold_id  <- basename(scaffold_dir)

cat("plotting SNP accuracy for:", scaffold_id, "\n")

# ---------- read accuracy-by-snp CSV ----------
csv_file <- file.path(scaffold_dir,paste0("accuracy_by_snp_", scaffold_id, ".csv"))

if (!file.exists(csv_file)) {stop("Missing accuracy_by_snp CSV for ", scaffold_id)}

accuracy_by_snp <- read.csv(csv_file, stringsAsFactors = FALSE)

# preserve SNP order as stored in CSV
accuracy_by_snp$SNP <- factor(accuracy_by_snp$SNP,levels = accuracy_by_snp$SNP)

# ---------- plot ----------
plot_snp_accuracy <- ggplot(accuracy_by_snp,
  aes(x = SNP, y = mean_abs_dif)) +
  geom_point(size = 0.6, alpha = 0.8) +
  theme_classic() +
  labs(title = paste0("SNP-wise mean absolute error: ", scaffold_id),
    x = "Position",
    y = "Mean absolute difference") +
  coord_cartesian( ylim = c(0, max(accuracy_by_snp$mean_abs_dif, na.rm = TRUE))) +
  theme(axis.text.x = element_blank())

# ---------- save ----------
out_png <- file.path(
  scaffold_dir,
  paste0("snp_accuracy_", scaffold_id, ".png"))

ggsave(filename = out_png,plot = plot_snp_accuracy,
  width = 8,
  height = 4,
  dpi = 300)

cat("Saved:", out_png, "\n")

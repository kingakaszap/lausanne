# 2026 jan-feb
# inspect accuracy of snipar imputed genotypes.

# libraries ----
.libPaths(c("/users/kkaszap/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(rhdf5)
library(SNPRelate) 
library(tidyverse) 

# step 1 - check scaffold folder
# if there is hdf5 file read it in
# -> if there is hdf5 file, also read in the corresponding vcf file
# the vcfs should already be filtered to only include the mums' id-s and genotypes.

# ---------- paths ----------
scaffolds_home <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/inferred_genotypes_by_scaffold"
vcf_home <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/scaffold_vcfs_filtered"
gds_home <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/per_scaffold_gds_files" 

# --------scaffold id extractor for matching ----
get_scaffold_id <- function(x) {sub("(_no_mum_chr_1.*$)|(_filtered.*$)", "", x)}

# ---------- get the slurm array id and find scaffold ----------
task_id <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
if (is.na(task_id)) stop("No SLURM_ARRAY_TASK_ID found")

scaffolds <- list.dirs(scaffolds_home, recursive = FALSE, full.names = TRUE)
scaffolds <- scaffolds[grepl("^Super-Scaffold_", basename(scaffolds))]

if (task_id > length(scaffolds)) {stop("Task ID exceeds number of scaffolds")}

scaffold_dir <- scaffolds[task_id]
folder_name <- basename(scaffold_dir)
scaffold_id <- get_scaffold_id(folder_name)

cat("Processing:", scaffold_id, "\n")

# ---- check if there is  HDF5 & read ----------
hdf_file <- list.files(scaffold_dir,pattern = "\\.hdf5$", full.names = TRUE)

if (length(hdf_file) == 0) {cat("No HDF5 for", scaffold_id, "\n")
  quit(save = "no")}

DosMat <- h5read(hdf_file, "imputed_par_gts")

pos    <- h5read(hdf_file, "pos")
ids    <- h5read(hdf_file, "families")

cat("inferred genotypes matrix' dimensions:", dim(DosMat), "\n")
cat("Length of ids:", length(ids), "\n")
cat("Length of pos:", length(pos), "\n")
cat("First 5 pos:", head(pos, 5), "\n")
cat("Last 5 pos:" , tail(pos, 5), "\n")
cat("First 5 ids:", head(ids, 5), "\n")

DosMat <- as.matrix(DosMat)
colnames(DosMat) <- as.character(ids)
rownames(DosMat) <- as.character(pos)

cat("HDF5 matrix:", dim(DosMat), "\n")
cat("First 5 SNPs in HDF5 matrix:",paste(head(rownames(DosMat), 5), collapse = " "), "\n")
cat("Last 5 SNPs in HDF5 matrix:",paste(tail(rownames(DosMat), 5), collapse = " "), "\n")

# ---------- read matching  GDS ----------

gds_file <- file.path(gds_home, paste0(scaffold_id, ".gds"))

if (!file.exists(gds_file)) {
  stop("GDS file missing for ", scaffold_id)}

gds <- snpgdsOpen(gds_file)
real_gts <- t(snpgdsGetGeno(gds))

sample_ids_real_gts  <- read.gdsn(index.gdsn(gds, "sample.id"))
snp_positions_real_gts <- read.gdsn(index.gdsn(gds, "snp.position"))
colnames(real_gts) <- as.character(sample_ids_real_gts)
rownames(real_gts) <- as.character(snp_positions_real_gts)

cat ("first 5 snps in vcf:" , head(snp_positions_real_gts, 5), "\n")
cat ("last 5 snps in vcf:" , tail(snp_positions_real_gts, 5), "\n")

cat("VCF matrix dimensions:", dim(real_gts), "\n")

snpgdsClose(gds)

# filter hdf5 to only include id-s that are also in the vcf, ie that we can compare -----
# individuals present in both
# (vcf has already been filtered in bash)
common_ids <- intersect(colnames(DosMat), sample_ids_real_gts)

cat("Individuals in DosMat (imputed parental gts):", ncol(DosMat), "\n")
cat("Individuals in VCF:", length(sample_ids_real_gts), "\n")
cat("Individuals present in both:", length(common_ids), "\n")
DosMat <- DosMat[, colnames(DosMat) %in% common_ids, drop = FALSE]
# keeps the duplicates

stopifnot(all(colnames(DosMat) %in% sample_ids_real_gts))

cat("Filtered DosMat dimensions:", dim(DosMat), "\n")
cat("unique owl id-s in dosmat:", length(unique(colnames(DosMat))), "\n")

stopifnot(identical(rownames(DosMat), rownames(real_gts))) # another sanity check

# -----accuracy ------
# get dosage difference matrix. Ideally keep the duplicates. Do row & col based matching
out_accuracy_dir<-"/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/accuracy_summaries"

# make out dirs
if (!dir.exists(out_accuracy_dir)) {
  dir.create(out_accuracy_dir, recursive = TRUE)}
scaffold_dir_out_accuracy <- file.path(out_accuracy_dir, scaffold_id)
if (!dir.exists(scaffold_dir_out_accuracy)) {
  dir.create(scaffold_dir_out_accuracy, recursive = TRUE)}

# make empty dif matrix
diff_matrix <- matrix(NA, nrow = nrow(DosMat), ncol = ncol(DosMat))
rownames(diff_matrix) <- as.character(rownames(DosMat))
colnames(diff_matrix) <- as.character(colnames(DosMat))

# Loop over columns of DosMat (duplicated individuals included)
# match SNPs that exist in both matrices - should be all
# might not have to make common snp-s
common_snps <- intersect(rownames(DosMat), rownames(real_gts))

# sanity check 
cat("\n--- SNP sanity: before diff_matrix ---\n")
cat("common_snps head:", head(common_snps, 10), "\n")
cat("common_snps tail:", tail(common_snps, 10), "\n")

stopifnot(identical(common_snps, rownames(DosMat)))
# basically common snp-s should be all snp-s always

# are SNP names aligned exactly in order? not sure if this warning is correct
if (!all(common_snps == rownames(real_gts)[match(common_snps, rownames(real_gts))])) {
  warning("SNPs in DosMat and real_gts are NOT in the same order! Check alignment.")}

# stop if wrong
stopifnot(all(common_snps == rownames(real_gts)[match(common_snps, rownames(real_gts))]))

for (j in seq_len(ncol(DosMat))) {
  sample_id <- colnames(DosMat)[j]
  # Only fill in if sample exists in real_gts - should be always
  if (sample_id %in% colnames(real_gts)) {

    diff_matrix[common_snps, j] <- DosMat[common_snps, j] - real_gts[common_snps, sample_id]} else {
    warning("Sample ", sample_id, " not found in real_gts, leaving NA") }}

cat("dif matrix:", dim(diff_matrix), "\n")

#  build the comparison dataframe safely
comparison_df <- data.frame(scaffold_id = scaffold_id,
  SNP         = rep(rownames(DosMat), times = ncol(DosMat)),
  Sample      = rep(colnames(DosMat), each = nrow(DosMat)),
  Predicted   = as.vector(DosMat),
  Actual      = as.vector(DosMat - diff_matrix),  # get the actual values
  abs_dif     = as.vector(abs(diff_matrix)),
  stringsAsFactors = FALSE)

# mean absolute difference
mean_abs_dif <- mean(comparison_df$abs_dif, na.rm = TRUE)

write.csv(comparison_df, file = file.path(scaffold_dir_out_accuracy, paste0("comparison_long_df_", scaffold_id, ".csv")),
          row.names = FALSE)

# comparison_df <- comparison_df_full %>% filter(!is.na(Predicted))
comparison_df$abs_dif_subtract <- abs(comparison_df$Predicted - comparison_df$Actual)

mean(comparison_df$abs_dif_subtract, na.rm = TRUE)

# sanity checks 
cat("\n--- SNP sanity: NA rows in diff_matrix ---\n")

na_rows <- which(rowSums(!is.na(diff_matrix)) == 0)

cat("Number of all-NA SNPs:", length(na_rows), "\n")

cat("First 10 all-NA SNPs:",
    rownames(diff_matrix)[na_rows][1:10], "\n")

cat("Last 10 all-NA SNPs:",
    tail(rownames(diff_matrix)[na_rows], 10), "\n")

# for comparison df

cat("\n--- SNP sanity: comparison_df ---\n")

cat("Unique SNPs in comparison_df:",
    length(unique(comparison_df$SNP)), "\n")

stopifnot(identical(
    sort(unique(comparison_df$SNP)),
    as.character(pos)))


# accuracy by snp-----


# Explicitly preserve HDF5 order and IDs
accuracy_by_snp <- comparison_df %>%
  # Keep all SNPs from DosMat in original HDF5 order
  mutate(SNP = factor(SNP, levels = rownames(DosMat))) %>%
  group_by(SNP) %>%
  summarise(
    mean_abs_dif = mean(abs_dif, na.rm = TRUE),
    n_non_na = sum(!is.na(abs_dif)),
    .groups = "drop"
  ) %>%
  arrange(SNP)
all(levels(accuracy_by_snp$SNP) == as.character(pos))


# Preserve SNP order from DosMat
accuracy_by_snp <- accuracy_by_snp %>%
  mutate(SNP = factor(SNP, levels = rownames(DosMat))) %>%
  arrange(SNP)
stopifnot(all(levels(accuracy_by_snp$SNP) == rownames(DosMat)))


cat("\n--- SNP sanity: accuracy_by_snp ---\n")

cat("Total SNPs:", nrow(accuracy_by_snp), "\n")
cat("NA mean_abs_dif:",
    sum(is.na(accuracy_by_snp$mean_abs_dif)), "\n")

cat("First 10 NA SNPs:",
    accuracy_by_snp$SNP[is.na(accuracy_by_snp$mean_abs_dif)][1:10], "\n")

# Save CSV for this scaffold
write.csv(accuracy_by_snp,
          file = file.path(scaffold_dir_out_accuracy, paste0("accuracy_by_snp_", scaffold_id, ".csv")),
          row.names = FALSE)

accuracy_by_snp <- accuracy_by_snp %>%arrange(SNP)

plot_snp_accuracy <- ggplot(accuracy_by_snp,
                            aes(x = SNP, y = mean_abs_dif)) +
  geom_point(size = 0.6, alpha = 0.8) +
  theme_classic() +
  labs(title = paste0("SNP-wise mean imputation error: ", scaffold_id),
    x = "SNP position",
    y = "Mean absolute difference") +
  coord_cartesian(ylim = c(0, max(accuracy_by_snp$mean_abs_dif, na.rm = TRUE))) +
  theme(axis.text.x = element_blank())

# individual wise mean accuracy -----
accuracy_by_ind <- comparison_df %>%
  group_by(Sample) %>%
  summarise(mean_abs_dif = mean(abs_dif, na.rm = TRUE),
    .groups = "drop")

# Save CSV
write.csv(accuracy_by_ind,
          file = file.path(scaffold_dir_out_accuracy, paste0("accuracy_by_individual_", scaffold_id, ".csv")),
          row.names = FALSE)


plot_ind_accuracy <- ggplot(accuracy_by_ind, aes(x = Sample, y = mean_abs_dif)) +
  geom_bar(stat = "identity") +
  theme_classic() +
  labs(title = paste0("Individual-wise mean accuracy: ", scaffold_id),
       x = "Sample",y = "Mean Absolute Difference") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
ggsave(filename = file.path(scaffold_dir_out_accuracy, paste0("individual_accuracy_", scaffold_id, ".png")),
       plot = plot_ind_accuracy,
       dpi = 300,width = 8,height = 4)

# accuracy avg per scaffold - 1row csv-s ----
scaffold_summary <- comparison_df %>%
  summarise(scaffold_id = first(scaffold_id),
            mean_abs_dif = mean(abs_dif, na.rm = TRUE),
            median_abs_dif = median(abs_dif, na.rm = TRUE),
            n_snps = n_distinct(SNP),n_samples = n_distinct(Sample))

# Save CSV
write.csv(scaffold_summary,
          file = file.path(out_accuracy_dir, paste0("scaffold_summary_", scaffold_id, ".csv")),
          row.names = FALSE)


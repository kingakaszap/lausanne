# 29/01/2026
# get NA metrics on all scaffolds from snipar imputation

# libraries ----
.libPaths(c("/users/kkaszap/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(rhdf5)
library(SNPRelate) 
library(tidyverse) 

# step 1 -
# check scaffold folder
# if there is hdf5 file read it in
# -> if there is hdf5 file, also read in the corresponding vcf file
# the vcfs should already be filtered to only include the mums' id-s and genotypes.
# ADD: filter the hdf5 to only include id-s that are present in the vcf/ gds.

# ---------- paths ----------
scaffolds_home <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/inferred_genotypes_by_scaffold"
vcf_home <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/scaffold_vcfs_filtered"
gds_home <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/per_scaffold_gds_files" 
# remember to create the out directory

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

# ---------- check if there is HDF5 ----------
hdf_file <- list.files(scaffold_dir,pattern = "\\.hdf5$", full.names = TRUE)

if (length(hdf_file) == 0) {cat("No HDF5 for", scaffold_id, "\n")
  quit(save = "no")}

# ---------- read HDF5 ----------
DosMat <- h5read(hdf_file, "imputed_par_gts")

pos    <- h5read(hdf_file, "pos")
ids    <- h5read(hdf_file, "families")

cat("DosMat dimensions:", dim(DosMat), "\n")
cat("Length of ids:", length(ids), "\n")
cat("Length of pos:", length(pos), "\n")
cat("First 5 pos:", head(pos, 5), "\n")
cat("First 5 ids:", head(ids, 5), "\n")

DosMat <- as.matrix(DosMat)
colnames(DosMat) <- ids
rownames(DosMat) <- as.character(pos)

cat("HDF5 matrix:", dim(DosMat), "\n")

# ---------- read matching VCF -> GDS ----------
vcf_file <- file.path(vcf_home,paste0(scaffold_id, "_filtered.vcf.gz"))

if (!file.exists(vcf_file)) {stop("no VCF for ", scaffold_id)}

gds_file <- file.path(gds_home, paste0(scaffold_id, ".gds"))

if (!dir.exists(gds_home)) {
  dir.create(gds_home, recursive = TRUE)}

if (!file.exists(gds_file)) {
    snpgdsVCF2GDS(vcf_file, gds_file, method = "biallelic.only", verbose = FALSE)}
gds <- snpgdsOpen(gds_file)
real_gts <- t(snpgdsGetGeno(gds))
colnames(real_gts) <- read.gdsn(index.gdsn(gds, "sample.id"))
rownames(real_gts) <- read.gdsn(index.gdsn(gds, "snp.position"))
sample_ids_real_gts <- read.gdsn(index.gdsn(gds, "sample.id"))
rsid <- read.gdsn(index.gdsn(gds, "snp.position"))
head(rsid, 5)
str(rsid)
  
cat("VCF matrix:", dim(real_gts), "\n")

snpgdsClose(gds)

# filter hdf5 to only include id-s that are also in the vcf, ie that we can compare -----
# individuals present in both
common_ids <- intersect(colnames(DosMat), sample_ids_real_gts)

cat("Individuals in DosMat (imputed parental gts):", ncol(DosMat), "\n")
cat("Individuals in VCF:", length(sample_ids_real_gts), "\n")
cat("Individuals present in both:", length(common_ids), "\n")
DosMat <- DosMat[, colnames(DosMat) %in% common_ids, drop = FALSE]

stopifnot(all(colnames(DosMat) %in% sample_ids_real_gts))

# quantify duplication
dup_counts <- table(colnames(DosMat))
summary(dup_counts)

cat("Filtered DosMat dimensions:", dim(DosMat), "\n")

# ---- % NA-s by individual and by SNP-s in scaffold ----
  # NA per ind
  na_summary_individual <-colMeans(is.na(DosMat))*100
  na_summary_individual_df <- data.frame(scaffold_id = scaffold_id,
                                         MumId = colnames(DosMat),
                              percent_missing_snps = na_summary_individual)
  
  # NP per SNP
  na_snp_summary<- rowMeans(is.na(DosMat))*100
  na_snp_df <- data.frame(position = rownames(DosMat),
                          scaffold_id = scaffold_id,
                          percent_missing_in_individuals = na_snp_summary)
  
out_na_dir<-"/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/na_summaries"
if (!dir.exists(out_na_dir)) {
  dir.create(out_na_dir, recursive = TRUE)}
# remember to make this dir
scaffold_dir_out <- file.path(out_na_dir, scaffold_id)
if (!dir.exists(scaffold_dir_out)) {
  dir.create(scaffold_dir_out, recursive = TRUE)}

saveRDS(na_summary_individual_df,
  file = file.path(scaffold_dir_out, paste0("na_summary_individual_df_", scaffold_id, ".rds")))
write.csv(na_summary_individual_df,file = file.path(scaffold_dir_out, paste0("na_summary_individual_df_", scaffold_id, ".csv")), row.names = FALSE)

saveRDS(na_snp_df,
  file = file.path(scaffold_dir_out, paste0("na_by_snp_", scaffold_id, ".rds")))
write.csv(na_snp_df,file = file.path(scaffold_dir_out, paste0("na_by_snp_", scaffold_id, ".csv")), row.names = FALSE)

plot_na_inds <- ggplot(na_summary_individual_df, aes(x = percent_missing_snps)) +
  geom_histogram(bins = 30, colour = "white") +
  theme_classic() +
  labs(title = scaffold_id,
    x = "Percentage of SNP-s missing",
    y = "Number of individuals\n") +
   theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))
ggsave(filename = file.path(scaffold_dir_out, paste0("missingness_byindividual_hist_", scaffold_id, ".png")),
  plot = plot_na_inds,dpi = 300, width = 6, height = 4)

plot_na_pos <-  ggplot(na_snp_df, aes(x = percent_missing_in_individuals)) +
  geom_histogram(bins = 30, colour = "white") +
  theme_classic() +
  labs(title = scaffold_id,
       x = "\n%of individuals in which SNP is missing",
       y = "Number of SNP-s\n") +
  theme(axis.title = element_text(size = 12),
        axis.text = element_text(size = 12))

ggsave(filename = file.path(scaffold_dir_out, paste0("missingness_by_snps_hist_", scaffold_id, ".png")),
       plot = plot_na_pos,dpi = 300, width = 6, height = 4)


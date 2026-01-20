
#libs ----
.libPaths(c("/users/kkaszap/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(rhdf5)

# scaffold folders ----
scaffolds_home <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/inferred_genotypes_by_scaffold"
scaffolds <- list.dirs(scaffolds_home, recursive = FALSE, full.names = TRUE)
scaffolds <- scaffolds[grepl("^Super-Scaffold_", basename(scaffolds))]

cat("Found scaffolds:\n")
print(basename(scaffolds))

# function to read one HDF5 scaffold ---

read_hdf5_snipar <- function(scaffold_folder){
  hdf_file <- list.files(scaffold_folder, pattern = "\\.hdf5$", ignore.case = TRUE, full.names = TRUE)
  
  if (length(hdf_file) == 0){
    warning("No HDF5 file found in folder: ", scaffold_folder)
    return(NULL)}
  
  cat("Reading file:", hdf_file, "\n")
  
  DosMat <- as.matrix(h5read(hdf_file, "imputed_par_gts"))
  pos    <- h5read(hdf_file, "pos")
  ids    <- h5read(hdf_file, "families")
  
  # sanity checks
  if (length(pos) != nrow(DosMat)){
    stop("Length of pos does not match number of rows in DosMat in file: ", hdf_file)}
  if (length(ids) != ncol(DosMat)){
    stop("Length of families does not match number of columns in DosMat in file: ", hdf_file)}
  
  colnames(DosMat) <- ids
  scaffold_id <- sub("(_no_mum.*$)", "", basename(scaffold_folder))
  rownames(DosMat) <- paste0(scaffold_id, ":", pos)
  
  return(DosMat)}

# test on 1 scaffold ----
test_scaffold <- scaffolds[8]
test_mat <- read_hdf5_snipar(test_scaffold)

cat("\nTest scaffold:", basename(test_scaffold), "\n")
print(dim(test_mat))
print(head(rownames(test_mat)))
print(head(colnames(test_mat)))
print(table(test_mat, useNA = "ifany"))

# read in all scaffolds ----

genotypes_list <- lapply(scaffolds, read_hdf5_snipar)
genotypes_list <- genotypes_list[!sapply(genotypes_list, is.null)]


# ensure all individuals match
stopifnot( all(sapply(genotypes_list, function(x) identical(colnames(x), colnames(genotypes_list[[1]])))))

# merge all scaffolds
DosMat_all <- do.call(rbind, genotypes_list)

# final checks ----
cat("\nFinal matrix dimensions:\n")
print(dim(DosMat_all))

cat("\nFirst 5 SNPs:\n")
print(head(rownames(DosMat_all), 5))

cat("\nFirst 5 individuals:\n")
print(head(colnames(DosMat_all), 5))

cat("\nDosage counts:\n")
print(table(DosMat_all, useNA = "ifany"))

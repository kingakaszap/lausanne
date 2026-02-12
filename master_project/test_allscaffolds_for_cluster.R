.libPaths(c("/users/kkaszap/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(rhdf5)

scaffolds_home <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/test_on_all_scaffolds/inferred_genotypes_by_scaffold"
scaffolds <- list.dirs(scaffolds_home, recursive = FALSE, full.names = TRUE)
# safety check
scaffolds <- scaffolds[grepl("^Super-Scaffold_", basename(scaffolds))]

# function to read each hdf5 ----
read_hdf5_snipar <- function(scaffolds){
  hdf_file <- list.files(scaffolds, 
                         pattern = "\\.hdf5$",
                         full.names = TRUE)
  if (length(hdf_file) == 0) {
    message("no hdf5 file in", scaffolds)
    return(NULL)}
  
  DosMat<-h5read(hdf_file, "imputed_par_gts")
  pos<- h5read(hdf_file, "pos")
  ids <- h5read(hdf_file, "families")
  
  DosMat <- as.matrix(DosMat)
  
  
  colnames(DosMat) <-ids
  scaffold_id <- sub("(_no_mum.*$)", "", basename(scaffolds))
  rownames(DosMat)<- paste0(scaffold_id, ":", pos)
  
  return(DosMat) }  

# test on 1 scaff -----
test_scaffold <- scaffolds[8]
test_mat <- read_hdf5_snipar(test_scaffold)

cat("Test scaffold:", basename(test_scaffold), "\n")
print(dim(test_mat))
print(head(rownames(test_mat)))
print(head(colnames(test_mat)))
print(table(test_mat, useNA = "ifany"))

# read in all scaffolds
genotypes_list <- lapply(scaffolds, read_hdf5_snipar)
genotypes_list <- genotypes_list[!sapply(genotypes_list, is.null)]

stopifnot(
  all(sapply(genotypes_list, function(x)
    identical(colnames(x), colnames(genotypes_list[[1]]))
  ))
)


DosMat_all <- do.call(rbind, genotypes_list)


# check final merged matrix ----
cat("Final matrix dimensions:\n")
print(dim(DosMat_all))

cat("\nFirst 5 SNPs:\n")
print(head(rownames(DosMat_all), 5))

cat("\nFirst 5 individuals:\n")
print(head(colnames(DosMat_all), 5))

cat("\nDosage counts:\n")
print(table(DosMat_all, useNA = "ifany"))

# new stuff -- run after checks ----
output_file <- file.path(scaffolds_home, "DosMat_all.rds")
saveRDS(DosMat_all, output_file)
cat("\nMerged matrix saved to:", output_file, "\n")

family_file <- file.path(scaffolds_home, "individual_ids_in_hdf5.txt")
writeLines(unique(colnames(DosMat_all)), family_file)
cat("individual IDs saved to:", family_file, "\n")

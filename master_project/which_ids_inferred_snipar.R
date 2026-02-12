# which id-s have been inferred with snipar?

# libraries ----
.libPaths(c("/users/kkaszap/R/x86_64-pc-linux-gnu-library/4.4", .libPaths()))
library(rhdf5)
library(SNPRelate) 
library(tidyverse) 

hdf_file <- "/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/genotype_inference_all_input_data/superscaffolds_inferred_genotypes/Super-Scaffold_1_chr_1/Super-Scaffold_1_no_mum_with_map.hdf5"

DosMat <- h5read(hdf_file, "imputed_par_gts")

pos    <- h5read(hdf_file, "pos")
ids    <- h5read(hdf_file, "families")

cat("inferred genotypes matrix' dimensions:", dim(DosMat), "\n")
cat("Length of ids:", length(ids), "\n")
cat("Length of pos:", length(pos), "\n")
cat("First 5 pos:", head(pos, 5), "\n")
cat("Last 5 pos:" , tail(pos, 5), "\n")
cat("First 5 ids:", head(ids, 5), "\n")

# write csv of id-s
# wrap ids into a data frame
ids_vec <- as.vector(ids)
ids_df <- tibble(ID = ids_vec) %>%count(ID, name = "n_occurrences")

# write to csv
write_csv(ids_df,"/work/FAC/FBM/DEE/jgoudet/barn_owl/kkaszap/snipar/genotype_inference_all_input_data/ids_inferred_snipar.csv")

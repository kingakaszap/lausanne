library(SNPRelate)
library(rhdf5)
library(gdsfmt)
library(knitr)
library(tidyverse)

h5ls("master_project/data/actual_gt_inference/filtered_superscaffold3.hdf5")

DosMat_genotype_inference = h5read("master_project/data/actual_gt_inference/filtered_superscaffold3.hdf5", "imputed_par_gts")
Pos_gti = h5read("master_project/data/actual_gt_inference/filtered_superscaffold3.hdf5", "pos")
Families_imputed_gti = h5read("master_project/data/actual_gt_inference/filtered_superscaffold3.hdf5", "families")
dim(DosMat_genotype_inference)
# name the cols and the rows
rownames(DosMat_genotype_inference) <- paste0(Pos_gti)
colnames(DosMat_genotype_inference) <- paste0(Families_imputed_gti)
View(DosMat_genotype_inference)

print(Families_imputed_gti)

na_summary_gti <-colMeans(is.na(DosMat_genotype_inference)*100)
na_summary_df_gti <- data.frame(MumId =names(na_summary_gti),
                            percent_missing_snps = na_summary_gti)
View(na_summary_df_gti)

id_to_plot <- "XF0041"  # put various id-s in here
# has many children w same father, hist looks good : M040764
# M043819 - 2 children, does not look good
# cant explicitly check who the imputation was based on though
col_id <- which(Families_imputed_gti == id_to_plot)
length(col_id)  

dosages_vector <- DosMat_genotype_inference[, col_id]
length(dosages_vector)

dosages_vector <- dosages_vector[!is.na(dosages_vector)]
length(dosages_vector) 

# hist of imputed values
hist(dosages_vector, breaks = 50, xlab = "Dosage" )

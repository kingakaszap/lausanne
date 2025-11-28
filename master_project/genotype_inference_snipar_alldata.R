library(SNPRelate)
library(rhdf5)
library(gdsfmt)
library(knitr)
library(tidyverse)

h5ls("master_project/data/actual_gt_inference/super_scaffold_3_all_data.hdf5")

DosMat_genotype_inference = h5read("master_project/data/actual_gt_inference/super_scaffold_3_all_data.hdf5", "imputed_par_gts")
Pos_gti = h5read("master_project/data/actual_gt_inference/super_scaffold_3_all_data.hdf5", "pos")
Families_imputed_gti = h5read("master_project/data/actual_gt_inference/super_scaffold_3_all_data.hdf5", "families")
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
length(unique(na_summary_df_gti$MumId)) # 221 individuals

id_to_plot <- "XF0086"  # put various id-s in here
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

# which snp-s are missing? 
head(DosMat_genotype_inference)

missingness_by_snp <- data.frame(
  snp_id = rownames(DosMat_genotype_inference),
  n_NA = rowSums(is.na(DosMat_genotype_inference)),
  percent_NA = rowMeans(is.na(DosMat_genotype_inference)) * 100
)

# View(missingness_by_snp)
# missing 
(ggplot (missingness_by_snp, aes(x = percent_NA))+
    geom_histogram(col = "white")+
    ggtitle("missingness by snp-s in inferred set")+
    theme_classic()+
    labs(y = "number of snps\n")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)))
ggsave("master_project/plots/snps_missing_inferenceset.png", dpi = 600)



sum(missingness_by_snp$n_NA>(ncol(DosMat_genotype_inference)*0.5)) # 4.7%
9712/nrow(DosMat_genotype_inference)
# 9712 out of ~ 2000 000 snp missing in over 50% of inds

sum(missingness_by_snp$n_NA<24)
# 79% of snp-s are present in >= 90% of individuals

# what cutoff???

# intersect with test set
# 1 - look at those missing in over 50% of individuals (worst ones)
missing_in_half<- missingness_by_snp %>% 
  filter(missingness_by_snp$n_NA>(ncol(DosMat_genotype_inference)*0.5))
nrow(missing_in_half) # ok

# relies on having made this df in other script !!!
missing_in_half_test<- missingness_by_snp_test %>% 
  filter(missingness_by_snp$n_NA>(ncol(DosMat_imputed)*0.5))
nrow(missing_in_half_test)

shared_values_bad_snps_test_and_gtinf <- intersect(missing_in_half$snp_id, missing_in_half_test$snp_id)
length(shared_values_bad_snps_test_and_gtinf)
# 8889 snp-s shared from 9712 --- basically all snp-s that are bad in test set 
# are bad in real gt inference; 
# and there are some extra new ones that are also bad :')

# good ones? the 80% present in > 90%
good_snps<- missingness_by_snp %>% 
  filter(n_NA<(ncol(DosMat_genotype_inference)*0.1))
nrow(good_snps)

good_snps_test <- missingness_by_snp_test %>% 
  filter(n_NA<(ncol(DosMat_imputed)*0.1))
nrow(good_snps_test)

shared_values_good_snps_test_and_gtinf <- intersect(good_snps$snp_id, good_snps_test$snp_id)
length(shared_values_good_snps_test_and_gtinf)
# 137161  - 83.7% shared
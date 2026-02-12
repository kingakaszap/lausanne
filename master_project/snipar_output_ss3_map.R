# libs ----
library(SNPRelate)
library(rhdf5)
library(gdsfmt)
library(knitr)
library(tidyverse)
# data ----
pedigree <- read_csv("master_project/data/coreccted_consensus_pedigree.txt")
sequenced <- read_csv("master_project/data/RingIds_seq1.txt")
sequenced_families_2_children <- read.csv("master_project/data/sequenced_families_2_children.csv")

# snipar output----
h5ls("master_project/data/Super-Scaffold_3_no_mum_with_map.hdf5")

DosMat_imputed_map = h5read("master_project/data/Super-Scaffold_3_no_mum_with_map.hdf5", "imputed_par_gts")
Pos_map = h5read("master_project/data/Super-Scaffold_3_no_mum_with_map.hdf5", "pos")
Families_imputed_map = h5read("master_project/data/Super-Scaffold_3_no_mum_with_map.hdf5", "families")

# name the cols and the rows
rownames(DosMat_imputed_map) <- paste0(Pos_map)
colnames(DosMat_imputed_map) <- paste0(Families_imputed_map)

print(Families_imputed_map) # mum id-s

counts <- c(zero = sum(DosMat_imputed_map == 0, na.rm = TRUE),
  one  = sum(DosMat_imputed_map == 1, na.rm = TRUE),
  two  = sum(DosMat_imputed_map == 2, na.rm = TRUE),
  miss = sum(!(DosMat_imputed_map %in% c(0,1,2))))
counts

head(DosMat_imputed_map)
# na check -----
na_summary_map <-colMeans(is.na(DosMat_imputed_map)*100)
na_summary_df_map <- data.frame(MumId =names(na_summary_map),
                            percent_missing_snps = na_summary_map)
# View(na_summary_df_map)

(ggplot(na_summary_df_map, aes(x = percent_missing_snps))+
    geom_histogram(col = "white")+
    ylab("Number of individuals\n")+
    xlab("\nMissing SNP-s (%)")+
    theme_classic()+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12)))
ggsave("master_project/plots/missing_snp-s_overall_with_map.png", dpi = 600)

children_number <- sequenced_families_2_children %>% group_by(MumId) %>% 
  mutate(in_families = n_distinct(DadId)) %>% 
  summarise(children = n(),
            in_families = first(in_families))
# View(children_number)

na_summary_df_map <- left_join(na_summary_df_map, children_number, by = "MumId")
# na_summary_df_filtered_map <- na_summary_df_map %>% 
#  filter(in_families == 1)

(scatter_na <- ggplot(na_summary_df_map, aes (x = children, y = percent_missing_snps))+
    geom_point()+
    theme_classic())
# accuracy -----
matrix_imputed_for_comparison_map<- DosMat_imputed_map[, colnames(DosMat_imputed_map)%in% sample_ids_real_gts]

# without removing dups
matrix_imputed_for_comparison_certain_map <- matrix_imputed_for_comparison_map
matrix_imputed_for_comparison_certain_map[!matrix_imputed_for_comparison_certain_map %in% c(0, 1, 2)] <- NA
# keeping only 1st duplicates as it didnt work with all of them
matrix_imputed_for_comparison_certain_map <- matrix_imputed_for_comparison_certain_map[, !duplicated(colnames(matrix_imputed_for_comparison_certain_map))]

common_samples <- intersect(colnames(matrix_imputed_for_comparison_certain_map),
                            colnames(gen_matrix_real_gts))
length(common_samples)  # how many samples actually have real data

# Create a vector of real genotypes named by SNP and Sample
y_names <- paste0(rep(rownames(gen_matrix_real_gts), times = ncol(gen_matrix_real_gts)),
  "_",rep(colnames(gen_matrix_real_gts), each = nrow(gen_matrix_real_gts)))
y_vector <- setNames(as.vector(gen_matrix_real_gts), y_names)

x_snp <- rep(rownames(matrix_imputed_for_comparison_certain_map), times = ncol(matrix_imputed_for_comparison_certain_map))
x_sample <- rep(colnames(matrix_imputed_for_comparison_certain_map), each = nrow(matrix_imputed_for_comparison_certain_map))

x_keys <- paste0(x_snp, "_", x_sample)
x_values <- as.vector(matrix_imputed_for_comparison_certain_map)

actual_values <- y_vector[x_keys]  # returns NA if the SNP+Sample is not in y_vector

comparison_df_certain_map <- data.frame(
  SNP = x_snp,
  Sample = x_sample,
  Predicted = x_values,
  Actual = actual_values,
  stringsAsFactors = FALSE)

comparison_df_certain_map$Match <- ifelse(is.na(comparison_df_certain_map$Predicted) | is.na(comparison_df_certain_map$Actual), NA,
                                          comparison_df_certain_map$Predicted == comparison_df_certain_map$Actual)
# View(comparison_df)
comparison_df_nona_map <- comparison_df_certain_map %>% filter(!is.na(Predicted))
# View(comparison_df_nona)

nrow(comparison_df_nona_map)
#42628491

sum(comparison_df_nona_map$Match == TRUE)
# 38869469
sum(comparison_df_nona_map$Match == FALSE)
# 3759022

38869469/42628491

# for cases when snipar is sure: 91% accuracy across all snp-s across all individuals. same

# Number of missing comparisons
num_missing <- sum(is.na(comparison_df_certain_map$Match))
num_missing # 25788429

conf_matrix <- comparison_df_nona_map %>% 
  count(Predicted, Actual) %>% 
  mutate(Predicted = as.factor(Predicted),
         Actual = as.factor(Actual)) %>% 
  pivot_wider(names_from=Actual, values_from=n,
              values_fill=0) %>% 
  arrange(Predicted)
colnames(conf_matrix) <- c("Predicted", "Actual_0", "Actual_1", "Actual_2")
conf_matrix        
conf_matrix <- conf_matrix %>% 
  mutate(Total = Actual_0 + Actual_2+ Actual_1)
kable(conf_matrix, caption = "Confusion Matrix: Predicted vs Actual Genotypes")


# percentage of incorrectly but surely imputed snp-s. 
summary_incorrect_map <- comparison_df_nona_map %>% 
  group_by(Sample) %>% 
  summarise(Percentage = sum(Match == TRUE)/n()*100)
# View(summary_incorrect)
summary_incorrect_map <- summary_incorrect_map %>% 
  rename(MumId = Sample)
View(summary_incorrect_map)
accuracy_families <- left_join(sequenced_families_2_children, summary_incorrect_map, by = "MumId")
View(accuracy_families)
write.table(accuracy_families, 
            "master_project/data/accuracy_families_with_map_ss3.csv",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
write.table(summary_incorrect_map, 
            "master_project/data/summary_incorrect_map_ss3_firstdupkept.csv",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
(histogram_incorrect <- ggplot(summary_incorrect_map, 
                               aes(x= Percentage))+
    geom_histogram(col = "white")+
    ggtitle("percent of correctly imputed snp-s from 'certain' values")+
    theme_classic()+
    labs(y = "number of individuals\n", x = "\nPercentage")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
    ))
ggsave("master_project/plots/correctness_only_certain_snps_map.png", dpi = 600)

# duplicates -----
matrix_imputed__with_dups_map<- DosMat_imputed_map[, colnames(DosMat_imputed_map)%in% sample_ids_real_gts]
dim(matrix_imputed__with_dups_map)

matrix_imputed_dups_map <- matrix_imputed__with_dups_map[, duplicated(colnames(matrix_imputed__with_dups_map)) 
                                                 | duplicated(colnames(matrix_imputed__with_dups_map), fromLast = TRUE)]
dim(matrix_imputed_dups_map)
# View(matrix_imputed_dups)

length(unique(colnames(matrix_imputed_dups_map)))
# 60 individuals are duplicated !!! ???
print(unique(colnames(matrix_imputed_dups_map)))

# NA values

na_summary_dups_map <-colMeans(is.na(matrix_imputed_dups_map)*100)
na_summary_dups_df_map <- data.frame(MumId =names(na_summary_dups_map),
                                 percent_missing_snps = na_summary_dups_map)
# View(na_summary_dups_df)
# up until now i think it is ok

# compare the different imputation accuracies 
# 1) avg difs
real_gts_dups <- gen_matrix_real_gts[, intersect(colnames(gen_matrix_real_gts), colnames(matrix_imputed_dups_map))]
length(unique(colnames(real_gts_dups))) # 60, all good

# make a dif matrix - but bear in mind that there are duplicates!
duplicates_all <- colnames(matrix_imputed_dups_map)
dim(matrix_imputed_dups_map)
base_names <- sub("\\.\\d+$", "", duplicates_all)
length(base_names)

# empty dif matrix
dif_matrix_dups_map <- 
  matrix(NA,
         nrow = nrow(matrix_imputed_dups_map),
         ncol = ncol(matrix_imputed_dups_map))
dim(dif_matrix_dups_map)
colnames(dif_matrix_dups_map) <- duplicates_all
rownames(dif_matrix_dups_map) <- rownames(matrix_imputed_dups_map)

for (i in seq_along(duplicates_all)) {
  if (base_names[i] %in% colnames(real_gts_dups)) {
    dif_matrix_dups_map[,i] <- real_gts_dups[, base_names[i]] - matrix_imputed_dups_map[,i]  
  } else {
    warning(paste("Column", base_names[i], "not found in real genotypes"))
  }
}

# View(dif_matrix_dups) 

colnames(matrix_imputed_dups_map) <- make.unique(colnames(matrix_imputed_dups_map))

real_matrix_expanded_dups_map <- real_gts_dups[, base_names, drop = FALSE]
head(real_matrix_expanded_dups_map)

stopifnot(all(dim(real_matrix_expanded_dups_map) == dim(matrix_imputed_dups_map)))

x_dups_map <- as.vector(matrix_imputed_dups_map)
y_dups_map <- as.vector(real_matrix_expanded_dups_map)

snp_labels_full_map <- rep(rownames(matrix_imputed_dups_map), times = ncol (matrix_imputed_dups_map))
individual_labels_full_map <- rep(colnames(matrix_imputed_dups_map), each = nrow(matrix_imputed_dups_map))

comparison_df_dups_map <- data.frame(
  SNP = snp_labels_full_map,
  MumId = individual_labels_full_map,
  Predicted = x_dups_map,
  Actual = y_dups_map,
  stringsAsFactors = FALSE
)

head(comparison_df_dups_map)

comparison_df_dups_map <- comparison_df_dups_map %>% 
  filter(!is.na(Predicted)) %>% 
  mutate(dif = Actual - Predicted)

summary_dups_map <- comparison_df_dups_map %>% 
  group_by(MumId) %>%  
  summarise (
    mean_abs_diff = mean(abs(dif))) %>% 
  ungroup()
View(summary_dups_map)

write.table(summary_dups_map, 
            "master_project/data/duplicates_with_map_ss3.csv",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)
# rounding accuracy DOES NOT WORK-----
# Create a vector of real genotypes named by SNP and Sample
y_names <- paste0(rep(rownames(gen_matrix_real_gts), times = ncol(gen_matrix_real_gts)),
                   "_",rep(colnames(gen_matrix_real_gts), each = nrow(gen_matrix_real_gts)))
y_vector <- setNames(as.vector(gen_matrix_real_gts), y_names)

x_snp_all <- rep(rownames(matrix_imputed_for_comparison_map), times = ncol(matrix_imputed_for_comparison_map))
x_sample_all <- rep(colnames(matrix_imputed_for_comparison_map), each = nrow(matrix_imputed_for_comparison_map))

x_keys_all <- paste0(x_snp_all, "_", x_sample_all)
x_values_all <- as.vector(matrix_imputed_for_comparison_map)

actual_values_all <- y_vector[x_keys_all]  # returns NA if the SNP+Sample is not in y_vector

comparison_df_map <- data.frame(
  SNP = x_snp_all,
  Sample = x_sample_all,
  Predicted = x_values_all,
  Actual = actual_values_all,
  stringsAsFactors = FALSE)



comparison_df_map <- comparison_df_map %>% 
  mutate(Rounded = round(Predicted))
comparison_df_map<- comparison_df_map %>% 
  mutate(match_rounded = (Actual == Rounded) )

summary_incorrect_rounded_map <- comparison_df_map %>% 
  group_by(Sample) %>% 
  summarise(Percentage = sum(match_rounded == TRUE)/n()*100)
 View(summary_incorrect_rounded_map)
summary_incorrect_rounded_map <- summary_incorrect_rounded_map %>% 
  rename(MumId = Sample)
(histogram_incorrect <- ggplot(summary_incorrect_rounded_map, 
                               aes(x= Percentage))+
    geom_histogram(col = "white")+
    ggtitle("percent of correctly imputed snp-s if values are rounded")+
    theme_classic()+
    labs(y = "number of individuals\n")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)))
ggsave("master_project/plots/correctness_rounded_values.png")
View(summary_incorrect_rounded)

# again, check relationship with number of children,
# this time with the rounded values
comparison_df_full<- comparison_df_full %>% 
  rename(MumId = Sample)

summary_incorrect_rounded <- left_join(summary_incorrect_rounded, children_number, by = "MumId")
summary_incorrect_rounded_filtered <- summary_incorrect_rounded %>% 
  filter(in_families == 1)
(scatter_rounded <- ggplot(summary_incorrect_rounded_filtered, aes (x = children, y = Percentage))+
    geom_point(colour = "purple")+
    labs(x = "\nNumber of chicks", y = "Percentage of SNP-s correct\n")+
    theme_classic()+
    theme(axis.text = element_text (size = 12),
          axis.title = element_text(size = 12)))
ggsave("master_project/plots/chicks_vs_correctness_with_rounded_values.png", dpi = 300)
cor(summary_incorrect_rounded_filtered$children, summary_incorrect_rounded_filtered$Percentage, method = "pearson", use = "complete.obs")


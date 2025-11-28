# inspect snipar output for SS3 for inferring mother genotypes 2025 sept-oct
# libs ----
# if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("SNPRelate")

library(SNPRelate)
library(rhdf5)
library(gdsfmt)
library(knitr)
library(tidyverse)


# to import ID-s that are in hdf5 file 
# DATA = snpgdsVCF2GDS("master_project/data/inferred_genotypes_from_snipar_with_id.vcf.gz", "master_project/data/inferred_genotypes_from_snipar_with_id.gds")

# REAL genotypes from vcf (to gds) of imputed mums-----
DATA = snpgdsOpen("master_project/data/inferred_genotypes_from_snipar_with_id.gds")
Samples = read.gdsn(index.gdsn(DATA, "sample.id"))

gen_matrix_real_gts <- snpgdsGetGeno('master_project/data/inferred_genotypes_from_snipar_with_id.gds')
dim(gen_matrix_real_gts)
colnames(gen_matrix_real_gts)
gen_matrix_real_gts<- t(gen_matrix_real_gts)
# View(gen_matrix_real_gts)
ls.gdsn(DATA)
head(read.gdsn(index.gdsn(DATA, "snp.position")), 10) # but these are the ones that correspond to the hdf5 positions.

snp_ids_real_gts<- read.gdsn(index.gdsn(DATA, "snp.rs.id")) # this is the snp id i need
sample_ids_real_gts<- read.gdsn(index.gdsn(DATA, "sample.id")) # and this is the id of individuals!
snp_pos_real_gts <- read.gdsn(index.gdsn(DATA,"snp.position" ))

length(sample_ids_real_gts)
dim(gen_matrix_real_gts)
rownames(gen_matrix_real_gts) <-paste0(snp_pos_real_gts)
colnames(gen_matrix_real_gts) <- paste0(sample_ids_real_gts)
# imputed genotypes of mums- output from snipar----

h5ls("master_project/data/Super-Scaffold_3_no_mum.hdf5")

DosMat_imputed = h5read("master_project/data/Super-Scaffold_3_no_mum.hdf5", "imputed_par_gts")
Pos = h5read("master_project/data/Super-Scaffold_3_no_mum.hdf5", "pos")
Families_imputed = h5read("master_project/data/Super-Scaffold_3_no_mum.hdf5", "families")

# name the cols and the rows
rownames(DosMat_imputed) <- paste0(Pos)
colnames(DosMat_imputed) <- paste0(Families_imputed)

print(Families_imputed) # mum id-s

counts <- c(
  zero = sum(DosMat_imputed == 0, na.rm = TRUE),
  one  = sum(DosMat_imputed == 1, na.rm = TRUE),
  two  = sum(DosMat_imputed == 2, na.rm = TRUE),
  miss = sum(!(DosMat_imputed %in% c(0,1,2)))
)
counts


# export family id-s (i.e. id-s of imputed moms) (should correspond to mum so that vcf can be filtered)
families_df <- data.frame(ID = Families_imputed)
length(unique(families_df$ID)) # 251 individuals imputed; ?
length(families_df$ID)
write.table ( families_df,
            file = "master_project/data/family_ids_from_hdf5.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)


# try to do histograms ???? -----

id_to_plot <- "M043819"  # put various id-s in here
# has many children w same father, hist looks good : M040764
# M043819 - 2 children, does not look good
# cant explicitly check who the imputation was based on though
col_id <- which(Families_imputed == id_to_plot)
length(col_id)  

dosages_vector <- DosMat_imputed[, col_id]
length(dosages_vector)

dosages_vector <- dosages_vector[!is.na(dosages_vector)]
length(dosages_vector) 

# hist of imputed values
hist(dosages_vector, breaks = 50, xlab = "Dosage" )
# families # for exchanging the initial id value


# check specific snp-s of duplicated individuals-----
sample_index <- which(sample_ids_real_gts == "M032352") # put any id in here
snp_index <- which(snp_ids_real_gts =="Super-Scaffold_3:6582188:A:C") # put snp of interest in here
 
real_genotype <- gen_matrix_real_gts[snp_index,sample_index] # correct dimensions, hopefully

real_genotype

# read in vcf of relatives ----
# i.e. individuals on whom imputation was based on.
# relatives_data = snpgdsVCF2GDS("master_project/data/Super-Scaffold_3_no_mum.vcf.gz", "master_project/data/Super-Scaffold_3_no_mum.gds")

relatives_data = snpgdsOpen("master_project/data/Super-Scaffold_3_no_mum.gds")
samples_relatives = read.gdsn(index.gdsn(relatives_data, "sample.id"))

gen_matrix_relatives <- snpgdsGetGeno('master_project/data/Super-Scaffold_3_no_mum.gds')
dim(gen_matrix_relatives)
gen_matrix_relatives <- t(gen_matrix_relatives)

ls.gdsn(relatives_data) # str-ish command
snp_ids_relatives<- read.gdsn(index.gdsn(relatives_data, "snp.rs.id")) # this is the snp id i need
sample_ids_relatives<- read.gdsn(index.gdsn(relatives_data, "sample.id"))
snp_ids_relatives_to_add_to_matrix <- read.gdsn(index.gdsn(relatives_data, "snp.position"))

rownames(gen_matrix_relatives) <- snp_ids_relatives_to_add_to_matrix
colnames(gen_matrix_relatives) <- sample_ids_relatives

sample_index <- which(sample_ids_relatives == "M032409") # put relative id here
snp_index <- which(snp_ids_relatives =="Super-Scaffold_3:1564204:C:A") # out snp of interest here

real_genotype <- gen_matrix_relatives[snp_index,sample_index] # the genotype. Real not necessary as there is no imputed genotype here
real_genotype 


# Basic statistics----
# 1.  % of NA values overall (hist & scatter as determined by no of children) ----

# imputed dataset: DosMat_imputed

# View(DosMat_imputed) # includes duplicates too.
na_summary <-colMeans(is.na(DosMat_imputed)*100)
na_summary_df <- data.frame(MumId =names(na_summary),
                            percent_missing_snps = na_summary)
# View(na_summary_df)

(ggplot(na_summary_df, aes(x = percent_missing_snps))+
  geom_histogram(col = "white")+
    ylab("Number of individuals\n")+
    xlab("\nMissing SNP-s (%)")+
    theme_classic()+
    theme(axis.title = element_text(size = 12),
          axis.text = element_text(size = 12)))
ggsave("master_project/plots/missing_snp-s_overall.png", dpi = 600)


# View(sequenced_families_2_children)
children_number <- sequenced_families_2_children %>% 
  group_by(MumId) %>% 
  mutate(in_families = n_distinct(DadId)) %>% 
  summarise(children = n(),
            in_families = first(in_families))
# View(children_number)

na_summary_df <- left_join(na_summary_df, children_number, by = "MumId")
na_summary_df_filtered <- na_summary_df %>% 
  filter(in_families == 1)

(scatter_na <- ggplot(na_summary_df_filtered, aes (x = children, y = percent_missing_snps))+
    geom_point()+
    theme_classic())
(scatter_na_all <- ggplot(na_summary_df, aes (x = children, y = percent_missing_snps))+
    geom_point()+
    theme_classic())

# 2.  % of incorrectly imputed snp-s-----
# - this will be harder. 
# first need to do the confusion matrix. 
# 1- filter imputed matrix to only include id-s that are 
# 1) present in the orignal gt-s matrix
# 2) -> from duplicately imputed, this will mean first imputation is included
# for now.

matrix_imputed_for_comparison<- DosMat_imputed[, colnames(DosMat_imputed)%in% sample_ids_real_gts]
# with duplicates --->

# back to only keeping 1 dup----
length(sample_ids_real_gts)
# remove dups except 1st occurence
matrix_imputed_for_comparison <- matrix_imputed_for_comparison[, !duplicated(colnames(matrix_imputed_for_comparison))]
dim(matrix_imputed_for_comparison) # should be 242 or less

# filter to snp-s where imputed values are 1,0 or 2
# replace uncertain to NA
matrix_imputed_for_comparison_certain <- matrix_imputed_for_comparison
matrix_imputed_for_comparison_certain[!matrix_imputed_for_comparison_certain %in% c(0, 1, 2)] <- NA

x <- as.vector(matrix_imputed_for_comparison_certain)
y <- as.vector(gen_matrix_real_gts)
snp_labels <- rep(rownames(matrix_imputed_for_comparison_certain), times = ncol(matrix_imputed_for_comparison_certain))
sample_labels <- rep(colnames(matrix_imputed_for_comparison_certain), each = nrow(matrix_imputed_for_comparison_certain))

comparison_df <- data.frame(
  SNP = snp_labels,
  Sample = sample_labels,
  Predicted = x,
  Actual = y,
  stringsAsFactors = FALSE
)

comparison_df$Match <- ifelse(is.na(comparison_df$Predicted) | is.na(comparison_df$Actual), NA,
                              comparison_df$Predicted == comparison_df$Actual)
# View(comparison_df)
comparison_df_nona <- comparison_df %>% filter(!is.na(Predicted))
# View(comparison_df_nona)

nrow(comparison_df_nona)
# 29 955 529

sum(comparison_df_nona$Match == TRUE)
# 27362731
sum(comparison_df_nona$Match == FALSE)
# 2592798

27362731/29955529

# for cases when snipar is sure: 91% accuracy across all snp-s across all individuals.

# Number of missing comparisons
num_missing <- sum(is.na(comparison_df$Match))
num_missing # 20 216 879

conf_matrix <- comparison_df_nona %>% 
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


comparison <- ifelse(is.na(x) | is.na(y), NA, paste0(x, "_", y))

conf_matrix <- table(comparison, useNA = "ifany")
conf_matrix

# percentage of incorrectly but surely imputed snp-s. 
summary_incorrect <- comparison_df_nona %>% 
  group_by(Sample) %>% 
  summarise(Percentage = sum(Match == TRUE)/n()*100)
# View(summary_incorrect)
summary_incorrect <- summary_incorrect %>% 
  rename(MumId = Sample)
(histogram_incorrect <- ggplot(summary_incorrect, 
                              aes(x= Percentage))+
  geom_histogram(col = "white")+
    ggtitle("percent of correctly imputed snp-s from 'certain' values")+
    theme_classic()+
    labs(y = "number of individuals\n", x = "\nPercentage")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)
          ))
ggsave("master_project/plots/correctness_only_certain_snps.png", dpi = 600)
# Correlation matrix -----
# start with matrix_imputed_for_comparison & gen_matrix_real_gts
diff_matrix <- matrix_imputed_for_comparison - gen_matrix_real_gts
# View(diff_matrix)

x_full <- as.vector(matrix_imputed_for_comparison)
y <- as.vector(gen_matrix_real_gts)
snp_labels_full <- rep(rownames(matrix_imputed_for_comparison), times = ncol(matrix_imputed_for_comparison))
sample_labels_full <- rep(colnames(matrix_imputed_for_comparison), each = nrow(matrix_imputed_for_comparison))

comparison_df_full <- data.frame(
  SNP = snp_labels_full,
  Sample = sample_labels_full,
  Predicted = x_full,
  Actual = y,
  stringsAsFactors = FALSE
)
comparison_df_full <- comparison_df_full %>% filter(!is.na(Predicted))
comparison_df_full$dif <- comparison_df_full$Predicted - comparison_df_full$Actual

mean(abs(comparison_df_full$dif))
# -0.09111837

comparison_df_full %>% 
  group_by(Actual) %>% 
  summarise(meandif = (mean(abs(dif))))


# check relationship between number of chicks used for imputation & imputation quality ------
# 1) extract number of children from dataset used to make the imputation dataset
# View(sequenced_families_2_children)
children_number <- sequenced_families_2_children %>% 
  group_by(MumId) %>% 
  mutate(in_families = n_distinct(DadId)) %>% 
  summarise(children = n(),
            in_families = first(in_families))
# View(children_number)
# View(summary_incorrect)
summary_incorrect
summary_incorrect <- left_join(summary_incorrect, children_number, by = "MumId")
summary_incorrect_for_plot  <- summary_incorrect %>% 
  filter(in_families == 1)
(scatter <- ggplot(summary_incorrect_for_plot, aes (x = children, y = Percentage))+
  geom_point(colour = "purple")+
    labs(x = "\nNumber of chicks", y = "Percentage of SNP-s correct\n")+
  theme_classic()+
    theme(axis.text = element_text (size = 12),
          axis.title = element_text(size = 12)))
ggsave("master_project/plots/chicks_vs_correctness_only_certain_values.png", dpi = 600)
# View(summary_incorrect_for_plot)
cor(summary_incorrect_for_plot$children, summary_incorrect_for_plot$Percentage, method = "pearson")


# with rounded values 



# hist with probabilities rounded to nearest even value -----
comparison_df_full <- comparison_df_full %>% 
  mutate(Rounded = round(Predicted))
comparison_df_full<- comparison_df_full %>% 
  mutate(match_rounded = (Actual == Rounded) )

summary_incorrect_rounded <- comparison_df_full %>% 
  group_by(Sample) %>% 
  summarise(Percentage = sum(match_rounded == TRUE)/n()*100)
# View(summary_incorrect_rounded)
summary_incorrect_rounded <- summary_incorrect_rounded %>% 
  rename(MumId = Sample)
(histogram_incorrect <- ggplot(summary_incorrect_rounded, 
                               aes(x= Percentage))+
    geom_histogram(col = "white")+
    ggtitle("percent of correctly imputed snp-s if values are rounded")+
    theme_classic()+
    labs(y = "number of individuals\n")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)))
ggsave("master_project/plots/correctness_rounded_values.png")

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


# look at duplicated individuals -----

matrix_imputed__with_dups<- DosMat_imputed[, colnames(DosMat_imputed)%in% sample_ids_real_gts]
dim(matrix_imputed__with_dups)

matrix_imputed_dups <- matrix_imputed__with_dups[, duplicated(colnames(matrix_imputed__with_dups)) 
                                                 | duplicated(colnames(matrix_imputed__with_dups), fromLast = TRUE)]
dim(matrix_imputed_dups)
# View(matrix_imputed_dups)

length(unique(colnames(matrix_imputed_dups)))
# 60 individuals are duplicated !!! ???
print(unique(colnames(matrix_imputed_dups)))

# NA values

na_summary_dups <-colMeans(is.na(matrix_imputed_dups)*100)
na_summary_dups_df <- data.frame(MumId =names(na_summary_dups),
                            percent_missing_snps = na_summary_dups)
# View(na_summary_dups_df)
# up until now i think it is ok

# compare the different imputation accuracies 
# 1) avg difs
real_gts_dups <- gen_matrix_real_gts[, intersect(colnames(gen_matrix_real_gts), colnames(matrix_imputed_dups))]
length(unique(colnames(real_gts_dups))) # 60, all good

# make a dif matrix - but bear in mind that there are duplicates!
duplicates_all <- colnames(matrix_imputed_dups)
dim(matrix_imputed_dups)
base_names <- sub("\\.\\d+$", "", duplicates_all)
length(base_names)

# empty dif matrix
dif_matrix_dups <- 
  matrix(NA,
         nrow = nrow(matrix_imputed_dups),
         ncol = ncol(matrix_imputed_dups))
dim(dif_matrix_dups)
colnames(dif_matrix_dups) <- duplicates_all
rownames(dif_matrix_dups) <- rownames(matrix_imputed_dups)

for (i in seq_along(duplicates_all)) {
  if (base_names[i] %in% colnames(real_gts_dups)) {
    dif_matrix_dups[,i] <- real_gts_dups[, base_names[i]] - matrix_imputed_dups[,i]  
  } else {
    warning(paste("Column", base_names[i], "not found in real genotypes"))
  }
}
  
# View(dif_matrix_dups) 

colnames(matrix_imputed_dups) <- make.unique(colnames(matrix_imputed_dups))

real_matrix_expanded_dups <- real_gts_dups[, base_names, drop = FALSE]
head(real_matrix_expanded_dups)

stopifnot(all(dim(real_matrix_expanded_dups) == dim(matrix_imputed_dups)))

x_dups <- as.vector(matrix_imputed_dups)
y_dups <- as.vector(real_matrix_expanded_dups)

snp_labels_full <- rep(rownames(matrix_imputed_dups), times = ncol (matrix_imputed_dups))
individual_labels_full <- rep(colnames(matrix_imputed_dups), each = nrow(matrix_imputed_dups))

comparison_df_dups <- data.frame(
  SNP = snp_labels_full,
  MumId = individual_labels_full,
  Predicted = x_dups,
  Actual = y_dups,
  stringsAsFactors = FALSE
)

head(comparison_df_dups)

comparison_df_dups <- comparison_df_dups %>% 
  filter(!is.na(Predicted)) %>% 
  mutate(dif = Actual - Predicted)

summary_dups <- comparison_df_dups %>% 
  group_by(MumId) %>%  
  summarise (
             mean_abs_diff = mean(abs(dif))) %>% 
  ungroup()
# View(summary_dups)

# try to link with na %
na_summary_unique_dups_df <- na_summary_dups_df
na_summary_unique_dups_df <- na_summary_dups_df %>%
  mutate(across(1, ~ make.unique(as.character(.x))))
# View(na_summary_unique_dups_df)
dups_na_accuracy<- left_join (summary_dups, na_summary_unique_dups_df, by = "MumId")
# View(dups_na_accuracy)
(misssingness_vs_accuracy_dups<-
    ggplot(dups_na_accuracy, aes(x = percent_missing_snps, y = mean_abs_diff )) + geom_point(color = "purple")+
    theme_classic()+
    labs (x= "\nPercentage of NA snp-s", y = "mean absolute difference\n")+
    theme(axis.text = element_text(size = 12),
          axis.title=element_text(size = 12)))
ggsave("master_project/plots/misssingness_vs_accuracy_dups.png", dpi = 600)

model <- lm(percent_missing_snps~mean_abs_diff, data =dups_na_accuracy )
summary(model)
plot(model)
cor(dups_na_accuracy$percent_missing_snps,dups_na_accuracy$mean_abs_diff )
# try to do an average for the dups at each snp need to double check ! -----
unique_base <- unique(base_names)
# pre-allocate matrix for averaged values
imputed_avg_mat <- matrix(NA,
                          nrow = nrow(matrix_imputed_dups),
                          ncol = length(unique_base))
rownames(imputed_avg_mat) <- rownames(matrix_imputed_dups)
colnames(imputed_avg_mat) <- unique_base

# loop over each base MumId and take row-wise mean across duplicates
for (b in unique_base) {
  cols <- which(base_names == b)
  imputed_avg_mat[, b] <- rowMeans(matrix_imputed_dups[, cols, drop = FALSE], na.rm = TRUE)
}
real_matrix_avg <- real_gts_dups[, unique_base, drop = FALSE]
stopifnot(all(dim(real_matrix_avg) == dim(imputed_avg_mat)))
x_avg <- as.vector(imputed_avg_mat)
y_avg <- as.vector(real_matrix_avg)

snp_labels_full <- rep(rownames(imputed_avg_mat), times = ncol(imputed_avg_mat))
individual_labels_full <- rep(colnames(imputed_avg_mat), each = nrow(imputed_avg_mat))

comparison_df_avg <- data.frame(
  SNP = snp_labels_full,
  MumId = individual_labels_full,
  Predicted = x_avg,
  Actual = y_avg,
  stringsAsFactors = FALSE
)

comparison_df_avg <- comparison_df_avg %>%
  filter(!is.na(Predicted)) %>%
  mutate(dif = Actual - Predicted)

# View(comparison_df_avg)
# add the avg dif as a column for duplicates ----
avg_base_names <- colnames(imputed_avg_mat)

# Compute difference matrix
dif_avg_mat <- real_matrix_avg - imputed_avg_mat

# Compute mean absolute difference per base MumId
mean_abs_dif_avg <- colMeans(abs(dif_avg_mat), na.rm = TRUE)  # named vector
names(mean_abs_dif_avg) <- avg_base_names

# Extract base MumIds from original summary table (duplicates included)
dup_base_names <- sub("\\.\\d+$", "", summary_dups$MumId)

# Add new column with averaged mean absolute difference
summary_dups$mean_abs_diff_avg <- mean_abs_dif_avg[dup_base_names]

# View(summary_dups)

# look at avg dif between mum and dad genotypes ----
# family to look at first: 
# M020884
# Dad 1: M026602
# Dad 2: M032415

# extract mum column from REAL matrix: 
M020884 <- real_gts_dups[, 'M020884']
Dad1_M026602 <- gen_matrix_relatives[,'M026602']
Dad2_M032415 <- gen_matrix_relatives[, 'M032415']

head(Dad1_M026602)

M020884_family <- data.frame(
  SNP_pos = rownames(gen_matrix_relatives),
  M020884 = M020884,
  Dad1_M026602 = Dad1_M026602,
  Dad2_M032415 = Dad2_M032415
)

# View(M020884_family)
M020884_family <- M020884_family %>% 
  mutate(dif_from_M026602 =abs(M020884- Dad1_M026602),
         dif_from_M032415 = abs(M020884-Dad2_M032415))

(avg_dif_M026602 <- mean(M020884_family$dif_from_M026602))
(avg_dif_M032415 <- mean(M020884_family$dif_from_M032415))
# :(((

# calculate n_heterozygous
for_heterozygosity_M020884 <- M020884_family[, c("M020884", "Dad1_M026602", "Dad2_M032415")]
(heterozygosity_sum_M020884 <- apply(for_heterozygosity_M020884, 2, function(x) prop.table(table(factor(x, levels = 0:2)))) )
(het_counts_M020884 <- apply(for_heterozygosity_M020884, 2, function(x) table(factor(x, levels = 0:2))))
# dad 1 
38867/(20619+38867+147838)  # 0.1874699
# dad 2
40942/(147254+40942+19128) # 0.197

# add offspring!  ------

M026757 <- gen_matrix_relatives[, 'M026757']
M026758 <- gen_matrix_relatives[, 'M026758']
M026759 <- gen_matrix_relatives[, 'M026759']
M026798 <- gen_matrix_relatives[, 'M026798']
M026799 <- gen_matrix_relatives[, 'M026799']

M020884_family_with_M026602 <- data.frame(
  SNP_pos = rownames(gen_matrix_relatives),
  M020884 = M020884,
  Dad1_M026602 = Dad1_M026602,
  M026757 = M026757,
  M026758 = M026758,
  M026759 = M026759,
  M026798 = M026798,
  M026799 = M026799
)

# MPAD
M020884_family_with_M026602_matrix <- M020884_family_with_M026602[,4:8]
M020884_family_with_M026602_matrix <- M020884_family_with_M026602_matrix %>% 
  select_if( is.numeric)
M020884_family_with_M026602_matrix<-  as.matrix(M020884_family_with_M026602_matrix)

# function to compute MPAD
mean_pairwise_abs_diff <- function(dos_mat, na.action = c("pairwise", "impute")) {
  na.action <- match.arg(na.action)
  if (na.action == "impute") {
    # simple mean imputation per SNP
    snp_means <- rowMeans(dos_mat, na.rm = TRUE)
    inds_na <- is.na(dos_mat)
    dos_mat[inds_na] <- rep(snp_means, times = ncol(dos_mat))[inds_na]
  }
  # get all column pairs
  n <- ncol(dos_mat)
  pairs <- combn(n, 2)
  pair_means <- numeric(ncol(pairs))
  for (k in seq_len(ncol(pairs))) {
    a <- pairs[1, k]; b <- pairs[2, k]
    # compute absolute diff per SNP, ignoring SNPs with NA in either sample
    diffs <- abs(dos_mat[, a] - dos_mat[, b])
    pair_means[k] <- mean(diffs, na.rm = TRUE)
  }
  mean(pair_means)   # average across pairs
}
# run (pairwise NA handling)
mpad <- mean_pairwise_abs_diff(M020884_family_with_M026602_matrix, na.action = "pairwise")
mpad
# 0.19
# proportion in [0,1]:
mpad_prop <- mpad / 2
mpad_prop

# other family 
M032403 <- gen_matrix_relatives[, 'M032403']
M032405 <- gen_matrix_relatives[, 'M032405']
M032406 <- gen_matrix_relatives[, 'M032406']
M032407 <- gen_matrix_relatives[, 'M032407']
M032408 <- gen_matrix_relatives[, 'M032408']
M020884_family_with_M032415 <- data.frame(
  SNP_pos = rownames(gen_matrix_relatives),
  M020884 = M020884,
  Dad2_M032415 = Dad2_M032415,
  M032403 = M032403,
  M032405 = M032405,
  M032406 = M032406,
  M032407 = M032407,
  M032408 = M032408
)

M020884_family_with_M03241_matrix <- M020884_family_with_M032415[,4:8]
M020884_family_with_M03241_matrix <- M020884_family_with_M03241_matrix %>% 
  select_if( is.numeric)
M020884_family_with_M03241_matrix<-  as.matrix(M020884_family_with_M03241_matrix)

mpad <- mean_pairwise_abs_diff(M020884_family_with_M03241_matrix, na.action = "pairwise")
mpad
# 0.209
# proportion in [0,1]:
mpad_prop <- mpad / 2
mpad_prop
  
# other dups to check ----
summary_dups <-summary_dups %>% 
  mutate(basename = sub("\\.\\d+$", "", MumId))
(summary_dups_to_check <- summary_dups %>% 
  group_by(basename) %>% 
  summarise(dif_range = max(mean_abs_diff)-min(mean_abs_diff)))
# View(summary_dups_to_check)

# try to see if rel between accuracy and missingness???----
str(summary_incorrect_rounded)
summary_incorrect_rounded <- left_join(summary_incorrect_rounded, na_summary_df, by = 'MumId')
(misssingness_vs_accuracy<-
  ggplot(summary_incorrect_rounded, aes(x = percent_missing_snps, y = Percentage )) + geom_point(color = "purple")+
    theme_classic()+
    labs (x= "\nPercentage of NA snp-s", y = "Percentage of correct snp-s (rounded)\n")+
    theme(axis.text = element_text(size = 12),
          axis.title=element_text(size = 12)))
ggsave("master_project/plots/misssingness_vs_accuracy.png", dpi = 600)
# WHY ? 
# missingness - is it similar to actually inferred individuals?----
missingness_by_snp_test <- data.frame(
  snp_id = rownames(DosMat_imputed),
  n_NA = rowSums(is.na(DosMat_imputed)),
  percent_NA = rowMeans(is.na(DosMat_imputed)) * 100
)
(ggplot (missingness_by_snp_test, aes(x = percent_NA))+
    geom_histogram(col = "white")+
    ggtitle("missingness by snp-s in test mums")+
    theme_classic()+
    labs(y = "number of snps\n")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)))
ggsave("master_project/plots/snps_missing_testmums.png", dpi = 600)

# View(missingness_by_snp_test)
# missing 

sum(missingness_by_snp_test$n_NA>(ncol(DosMat_imputed)*0.5)) # 5%
# 10953 out of 2000 000 snp missing in over 50% of inds
8890/nrow(DosMat_imputed) # 4.2 of snp-s missing in over 50% inds
sum(missingness_by_snp_test$n_NA<(ncol(DosMat_imputed)*0.1))
# 78% of snp-s are present in >= 90% of individuals
161332/nrow(DosMat_imputed) 
# what cutoff???
# # try everything with removing all duplicates ------
dup_cols <- names(which(table(colnames(DosMat_imputed)) > 1))
length(dup_cols)
# keep only columns that are NOT duplicated
matrix_no_dups <- DosMat_imputed[, !colnames(DosMat_imputed) %in% dup_cols]
# and have gt-s
# Find the intersection
common_samples <- intersect(colnames(matrix_no_dups), colnames(gen_matrix_real_gts))
length(common_samples)  # how many match
matrix_no_dups_sub <- matrix_no_dups[, common_samples]
gen_matrix_real_gts_sub <- gen_matrix_real_gts[, common_samples]

# flatten and make comparison
x_nodups <- as.vector(matrix_no_dups_sub)
y_nodups <- as.vector(gen_matrix_real_gts_sub)

snp_labels_nodups <- rep(rownames(matrix_no_dups_sub), times = ncol(matrix_no_dups_sub))
sample_labels_nodups <- rep(colnames(matrix_no_dups_sub), each = nrow(matrix_no_dups_sub))

comparison_df_nodups <- data.frame(
  SNP = snp_labels_nodups,
  Sample = sample_labels_nodups,
  Predicted = x_nodups,
  Actual = y_nodups,
  stringsAsFactors = FALSE
)
comparison_df_nodups <- comparison_df_nodups %>% filter(!is.na(Predicted))
comparison_df_nodups$dif <- comparison_df_nodups$Predicted - comparison_df_nodups$Actual

comparison_df_nodups <- comparison_df_nodups %>% 
  mutate(Rounded = round(Predicted)) %>% 
  mutate(match_rounded =(Actual == Rounded))

summary_incorrect_rounded_nodups <- comparison_df_nodups %>% 
  group_by(Sample) %>% 
  summarise(Percentage = sum(match_rounded == TRUE)/n()*100)
# View(summary_incorrect_rounded)
summary_incorrect_rounded_nodups <- summary_incorrect_rounded_nodups %>% 
  rename(MumId = Sample)
(histogram_incorrect_nodups <- ggplot(summary_incorrect_rounded_nodups, 
                                      aes(x= Percentage))+
    geom_histogram(col = "white")+
    ggtitle("percent of correctly imputed snp-s if values are rounded")+
    theme_classic()+
    labs(y = "number of individuals\n")+
    theme(plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 12),
          axis.title = element_text(size = 12)))
ggsave("master_project/plots/correctness_rounded_values_nodups.png")

# thank fuck

summary_incorrect_rounded_nodups <- left_join(summary_incorrect_rounded_nodups, children_number, by = "MumId")
head(summary_incorrect_rounded_nodups)
(scatter <- ggplot(summary_incorrect_rounded_nodups, aes (x = children.x, y = Percentage))+
    geom_point(colour = "purple")+
    labs(x = "\nNumber of chicks", y = "Percentage of SNP-s correct\n")+
    theme_classic()+
    theme(axis.text = element_text (size = 12),
          axis.title = element_text(size = 12)))
# fab
ggsave("master_project/plots/chicks_vs_correctness_only_certain_values_nodups.png", dpi = 600)

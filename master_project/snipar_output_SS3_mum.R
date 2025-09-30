# inspect snipar output for SS3 for inferring mother genotypes

# if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("SNPRelate")

library(SNPRelate)

# need to change !! to import ID-s that are in hdf5 file
DATA = snpgdsVCF2GDS("master_project/data/inferred_genotypes_from_snipar_with_id.vcf.gz", "master_project/data/inferred_genotypes_from_snipar_with_id.gds")

# change from this line
DATA = snpgdsOpen("master_project/data/inferred_genotypes_from_snipar_with_id.gds")
Samples = read.gdsn(index.gdsn(DATA, "sample.id"))

gen_matrix <- snpgdsGetGeno('master_project/data/inferred_genotypes_from_snipar_with_id.gds')
dim(gen_matrix)

ls.gdsn(DATA)
snp_ids<- read.gdsn(index.gdsn(DATA, "snp.rs.id")) # this is the snp id i need
sample_ids<- read.gdsn(index.gdsn(DATA, "sample.id")) # and this is the id of individuals!

# BiocManager::install("rhdf5")

library(rhdf5)
h5ls("master_project/data/Super-Scaffold_3_no_mum.hdf5")

DosMat = h5read("master_project/data/Super-Scaffold_3_no_mum.hdf5", "imputed_par_gts")
Pos = h5read("master_project/data/Super-Scaffold_3_no_mum.hdf5", "pos")
Families = h5read("master_project/data/Super-Scaffold_3_no_mum.hdf5", "families")

View(Families)
print(Families)

counts <- c(
  zero = sum(DosMat == 0, na.rm = TRUE),
  one  = sum(DosMat == 1, na.rm = TRUE),
  two  = sum(DosMat == 2, na.rm = TRUE),
  miss = sum(!(DosMat %in% c(0,1,2)))
)
counts




# export family id-s (should correspond to mum so that vcf can be filtered)
families_df <- data.frame(ID = Families)
length(unique(families_df$ID)) # 251 ??????



write.table ( families_df,
            file = "master_project/data/family_ids_from_hdf5.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)


# try to do histograms ????


id_to_plot <- "M043819" # has many children w same father, hist looks good : M040764
# M043819 - 2 children, does not look good
# cant explicitly check who the imputation was based on though
col_id <- which(Families == id_to_plot)
length(col_id)  

dosages_vector <- DosMat[, col_id]
length(dosages_vector)

dosages_vector <- dosages_vector[!is.na(dosages_vector)]
length(dosages_vector) 


hist(dosages_vector, breaks = 50, xlab = "Dosage" )


families


# check specific snp-s of duplicated individuals
sample_index <- which(sample_ids == "M032352")
snp_index <- which(snp_ids =="Super-Scaffold_3:6582188:A:C")

real_genotype <- gen_matrix[sample_index,snp_index]
real_genotype

# read in vcf of relatives

relatives_data = snpgdsVCF2GDS("master_project/data/Super-Scaffold_3_no_mum.vcf.gz", "master_project/data/Super-Scaffold_3_no_mum.gds")

relatives_data = snpgdsOpen("master_project/data/Super-Scaffold_3_no_mum.gds")
samples_relatives = read.gdsn(index.gdsn(relatives_data, "sample.id"))

gen_matrix_relatives <- snpgdsGetGeno('master_project/data/Super-Scaffold_3_no_mum.gds')
dim(gen_matrix_relatives)

ls.gdsn(relatives_data)
snp_ids_relatives<- read.gdsn(index.gdsn(relatives_data, "snp.rs.id")) # this is the snp id i need
sample_ids_relatives<- read.gdsn(index.gdsn(relatives_data, "sample.id")) 


sample_index <- which(sample_ids_relatives == "M032409")
snp_index <- which(snp_ids_relatives =="Super-Scaffold_3:1564204:C:A")

real_genotype <- gen_matrix_relatives[sample_index,snp_index]
real_genotype

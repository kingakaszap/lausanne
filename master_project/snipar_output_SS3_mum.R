# inspect snipar output for SS3 for inferring mother genotypes

# if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# BiocManager::install("SNPRelate")

library(SNPRelate)

# need to change !! to import ID-s that are in hdf5 file
DATA = snpgdsVCF2GDS("master_project/data/Super-Scaffold_3_no_mum.vcf.gz", "master_project/data/Super-Scaffold_3_no_mum.gds")
DATA = snpgdsOpen("Super-Scaffold_3_no_mum.gds")
Samples = read.gdsn(index.gdsn(DATA, "sample.id"))

gen_matrix <- snpgdsGetGeno('Super-Scaffold_3_no_mum.gds')
dim(gen_matrix)




# BiocManager::install("rhdf5")

library(rhdf5)

DosMat = h5read("master_project/data/Super-Scaffold_3_no_mum.hdf5", "imputed_par_gts")
Pos = h5read("~/Downloads/Super-Scaffold_3_no_mum.hdf5", "pos")
Families = h5read("~/Downloads/Super-Scaffold_3_no_mum.hdf5", "families")
# relatedness of parents 
# for expected inbreeding coefficient 

library(tidyverse)
library(hierfstat)
library(SNPRelate)
library(rhdf5)
library(gdsfmt)


parents_df <- read.csv("master_project/data/both_parents_not_all_chicks.csv")
pedigree<-read.table("master_project/data/coreccted_consensus_pedigree.txt", sep= "")

# extract and save mum and dad ids so that vcf can be filtered
parent_ids_for_grm <- unique(c(parents_df$MumId,parents_df$DadId))
length(parent_ids_for_grm) # 561

write.table(parent_ids_for_grm,file = "master_project/data/parent_ids_for_grm.txt",
            quote = FALSE,row.names = FALSE,col.names = FALSE)

# DATA = snpgdsVCF2GDS("master_project/data/filtered_vcf_grm_parents.vcf.gz", "master_project/data/parent_gts_for_grm.gds")
 
DATA = snpgdsOpen("master_project/data/parent_gts_for_grm.gds")
samples = read.gdsn(index.gdsn(DATA, "sample.id"))
rownames(gen_matrix_parents) <- samples
length (samples) #561 ok
gen_matrix_parents <- snpgdsGetGeno('master_project/data/parent_gts_for_grm.gds')
dim(gen_matrix_parents)
ls.gdsn(DATA)


# grm

parents_in_gds <- intersect(parent_ids_for_grm, samples)
length(parents_in_gds)  # great

# subset GDS to parent samples
filtered_matrix_parents  <- gen_matrix_parents[parents_in_gds, , drop = FALSE]




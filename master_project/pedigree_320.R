# find families that have 2 or more genotyped children and both parents genotyped
# after -> find grandparents etc like chains

library(tidyverse)
library(readr)

pedigree <- read_csv("master_project/data/coreccted_consensus_pedigree.txt")
#View(pedigree)
nrow(pedigree) # 11526
str(pedigree)
# add FamilyID column to pedigree
pedigree$FamilyID <- paste(pedigree$MumId, pedigree$DadId, sep = "_")

sequenced <- read_csv("master_project/data/RingIds_seq1.txt")
# View(sequenced)
str(sequenced)
nrow(sequenced) # 3085

# compare sequenced with vcf - should correspond
ids_in_vcf <- readLines("master_project/data/ids_in_vcf.txt")
ids_in_vcf <- data.frame( ID = ids_in_vcf)
#View(ids_in_vcf)
nrow(ids_in_vcf)

# check correspondance
all(sequenced$RingId %in% ids_in_vcf$ID) # true so ... ? 

# families with >= 2 sequenced children -----
full_families <- pedigree %>%
  filter(!is.na(MumId) & !is.na(DadId))
nrow(full_families) # 8506

# sequenced parents
sequenced_parents <- full_families %>% 
  filter(MumId %in% sequenced$RingId & DadId %in% sequenced$RingId)
nrow(sequenced_parents)
# 3830

# families of 3, where everyone is sequenced
sequenced_child <- sequenced_parents %>% 
  filter(RingId %in% sequenced$RingId)
nrow(sequenced_child)
# 1767
#View(sequenced_child)
# safety check - all they are in the vcf?
all_ids_seq_child <- unique(c(sequenced_child$RingId, 
                    sequenced_child$MumId, 
                    sequenced_child$DadId))

all(all_ids_seq_child %in% ids_in_vcf$ID) # true SO HOW does it become false later

# how many children associated with each sequenced parent from these families?
mum_counts <- table(sequenced_child$MumId) 
dad_counts <- table(sequenced_child$DadId)

multiple_mums <- names(mum_counts[mum_counts >= 2])
multiple_dads <- names(dad_counts[dad_counts >= 2])

sequenced_families_2_children <- sequenced_child %>%
  filter(MumId %in% multiple_mums & DadId %in% multiple_dads)
# potential issue: what if they appear multiple times
# but with different partners?

nrow(sequenced_families_2_children)
# 1701
#View(sequenced_families_2_children)
# safety check again 
all_ids_multiple_children <- unique(c(sequenced_families_2_children$RingId, 
                                      sequenced_families_2_children$MumId, 
                                      sequenced_families_2_children$DadId))
all (all_ids_multiple_children %in% ids_in_vcf$ID) # still true


(summary <- sequenced_families_2_children %>% 
  group_by(MumId, DadId) %>% 
    tally(name = "Children"))
#View(summary)
# 410

length(unique(sequenced_families_2_children$MumId)) # 303

test_data_mom_removed <- sequenced_families_2_children 
test_data_mom_removed$MumId <- NA 
  
# text file of ID-s - with mother removed (remember that in pedigree, she needs to be kept)
ids_test_data_no_mum <- unique(unlist(test_data_mom_removed))
ids_test_data_no_mum <- ids_test_data_no_mum[!is.na(ids_test_data_no_mum)]
nomum_ids <- data.frame( ID = ids_test_data_no_mum )
#View(nomum_ids)
# check again with vcf 
all(nomum_ids$ID%in% ids_in_vcf$ID) 
# as per, issues were the NA-s


# write df
write.table(nomum_ids, 
            "master_project/data/test_ids_no_mum.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

# triple check
final_nomum_ids <- read.table("master_project/data/test_ids_no_mum.txt",
                             header = FALSE,
                             stringsAsFactors = FALSE)                         

#View(final_nomum_ids)
all(final_nomum_ids$V1%in% ids_in_vcf$ID) 
# amazing i can finally move on :')

# just mom id-s so i can have a separate vcf
test_data_just_mom_ids <- sequenced_families_2_children 
test_data_just_mom_ids$RingId <- NA 
test_data_just_mom_ids$DadId <- NA
test_data_just_mom_ids$FamilyID <- NA
#View(test_data_just_mom_ids)

test_data_just_mom_ids <- unique(unlist(test_data_just_mom_ids))
test_data_just_mom_ids <- test_data_just_mom_ids[!is.na(test_data_just_mom_ids)]
mum_ids <- data.frame( ID = test_data_just_mom_ids )
#View(mum_ids)
# check again with vcf 
all(mum_ids$ID%in% ids_in_vcf$ID) 
# great great
# but actually i would need to extract all the id-s from the hdf5 file !

# write df
write.table(mum_ids, 
            "master_project/data/test_mum_ids_intended.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)



# 2 children 1 parent? ----
# mum sequenced, 2 children or more also sequenced, dad not sequenced----
mums <- pedigree %>% 
  filter(MumId %in% sequenced$RingId) %>% 
  filter(RingId %in% sequenced$RingId) %>% 
  filter(!DadId %in% sequenced$RingId ) %>% 
  na.omit()

mum_counts2 <- table(mums$MumId)
# View(mum_counts2)

multiple_mums_2 <- names(mum_counts2[mum_counts2 >= 2])
length(multiple_mums_2) # 71

filtered_mums <- mums %>%
  filter(MumId %in% multiple_mums_2)
# View(filtered_mums) 

# how many families with this setup

summary_mums <- filtered_mums %>%
  group_by(MumId, DadId) %>%
  tally(name = "count")
# View(summary_mums) # 77-2 = 75 why is it not the same as 71 before whoooo

# dad sequenced, chicks too, mom no ----
dads <- pedigree %>% 
  filter(DadId %in% sequenced$RingId) %>% 
  filter(RingId %in% sequenced$RingId) %>% 
  filter(!MumId %in% sequenced$RingId) %>% 
  na.omit()

dad_counts2 <- table(dads$DadId)
multiple_dads2 <- names(dad_counts2[dad_counts2 >= 2])
filtered_dads <- dads %>% 
  filter(DadId %in% multiple_dads2)
# View(filtered_dads)

summary_dads <- filtered_dads %>%
  group_by(MumId, DadId) %>%
  tally(name = "count")
# View(summary_dads) # 56-7 = 49 


# where both parents missing, but 2+ children? -----
no_parents <- pedigree %>% 
  filter(!DadId %in% sequenced$RingId) %>% 
  filter(RingId %in% sequenced$RingId) %>% 
  filter(!MumId %in% sequenced$RingId) %>% 
  na.omit()
# View(no_parents)
summary_noparents <- no_parents %>% 
  group_by(FamilyID) %>% 
  tally(name = "count") %>% 
  filter(count >1)
no_parents<- no_parents %>% filter(FamilyID %in% summary_noparents$FamilyID)

# View(summary_noparents) # 47

# how many unsequenced children in the pedigree? -----
unsequenced_children <- pedigree %>% 
  filter(!RingId %in% sequenced$RingId)
nrow(unsequenced_children)
nrow(unsequenced_children)/nrow(pedigree) # wow

sequenced_children <- pedigree %>% 
  filter(RingId %in% sequenced$RingId)
nrow(sequenced_children) # 3085


# look at families where both parents, but not all siblings ----

# sequenced IDs as a character vector
sequenced_ids <- as.character(sequenced$RingId)

# summarise per family
family_summary <- pedigree %>%
  group_by(FamilyID, MumId, DadId) %>% # just for safety - should be the same grouping
  summarise(
    num_children = n(),
    num_sequenced_children = sum(RingId %in% sequenced_ids),
    num_unsequenced_children = sum(!(RingId %in% sequenced_ids)),
    both_parents_sequenced = all(c(MumId, DadId) %in% sequenced_ids),
    .groups = "drop"
  )
#View(family_summary)

not_all_children <- family_summary %>%
  filter(
    both_parents_sequenced,
    num_children > 1,
    num_sequenced_children > 0,
    num_unsequenced_children > 0
  )
nrow(not_all_children) # 354

# get relatives of to-be-imputed individuals ----
inds_to_impute <- rbind(filtered_mums, filtered_dads,no_parents )
# View(inds_to_impute)
all(inds_to_impute$RingId %in% c(filtered_mums$RingId, filtered_dads$RingId, no_parents$RingId))
# i hope all good
inds_to_impute <- inds_to_impute %>%  select(-FamilyID)
# extract individuals ID-s
inds_to_impute_vector <- unique(unlist(inds_to_impute, use.names = FALSE))
inds_to_impute_df <- data.frame(ID = inds_to_impute_vector)
# View(inds_to_impute_df)
write.table(inds_to_impute_df, 
            "master_project/data/ids_for_imputation_all.txt",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

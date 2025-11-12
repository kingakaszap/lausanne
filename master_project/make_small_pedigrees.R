library(readr)
pedigree_snipar <- read.delim("master_project/data/snipar_pedigree.txt", sep = " ")
View(pedigree_snipar)
pedigree_M020884 <- filter(pedigree_snipar, MOTHER_ID == "M020884")
View(pedigree_M020884)


# family 1
pedigree_M020884_with_M026602 <- filter(pedigree_M020884, FATHER_ID == "M026602")
View(pedigree_M020884_with_M026602)

write.table(pedigree_M020884_with_M026602, 
            file = "master_project/data/test_duplicates_M020884/pedigree_M020884_with_M026602.txt", 
            sep = " ",          
            row.names = FALSE,  
            col.names = TRUE,   
            quote = FALSE) 


# family 2
pedigree_M020884_with_M032415 <- filter(pedigree_M020884, FATHER_ID == "M032415")
pedigree_M020884_with_M032415 <- filter(pedigree_M020884_with_M032415, IID != "M032404")
View(pedigree_M020884_with_M032415)

write.table(pedigree_M020884_with_M032415, 
            file = "master_project/data/test_duplicates_M020884/pedigree_M020884_with_M032415.txt", 
            sep = " ",          
            row.names = FALSE,  
            col.names = TRUE,   
            quote = FALSE) 

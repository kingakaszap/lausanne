library(tidyverse)
refpanel <- read.table("master_project/data/REFPANEL_snps_LG_coords.txt")

str(refpanel)
refpanel <- refpanel %>%  
  rename(chrom = V1,
         pos = V2,
         LG = V3,
         newpos = V4)

scaffold3 <- filter(refpanel, chrom == "Super-Scaffold_3" )
str(scaffold3)
range(scaffold3$newpos)
range(scaffold3$pos)

scaffold49<- filter(refpanel, chrom == "Super-Scaffold_49")
str(scaffold49)

range(scaffold49$newpos)
range(scaffold49$pos)

length(unique(refpanel$chrom))


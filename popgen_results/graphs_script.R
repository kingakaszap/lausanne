library(tidyverse)

# base folder containing the 8 subfolders
base_path <- "popgen_results/Forward_Slim_migrants_from_GSE/"

# list all csv files inside all subfolders
all_files <- list.files(base_path, 
                        pattern = "\\.csv$", 
                        recursive = TRUE, 
                        full.names = TRUE)

# read all files
all_data <- map_df(all_files, function(f) {
  
  # extract file name only
  fname <- basename(f)
  
  # extract M, K, r using regex
  M <- str_extract(fname, "M=\\d+") %>% str_remove("M=") %>% as.numeric()
  K <- str_extract(fname, "KF=\\d+") %>% str_remove("KF=") %>% as.numeric()
  r <- str_extract(fname, "r\\d+") %>% str_remove("r") %>% as.numeric()
  
  # read the file
  df <- read_csv(f)
  
  # attach identifiers
  df %>%
    mutate(M = M, K = K, r = r)
})

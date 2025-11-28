library(tidyverse)

# base folder containing the 8 subfolders
base_path <- "./popgen_results/Forward_Sim_migrants_from_GSE/"

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
head(all_data)
str(all_data)
View(all_data)
str(all_data)
unique(all_data$M)
all_data<-data.frame((all_data))
all_data <- all_data %>%
  mutate(Year = map_int(Year, as.integer))

plot_data_froh_l100 <- all_data %>%
  filter(K == 100,
         Year == 55702)
unique(plot_data_froh_l100$M)

(plot <- ggplot(plot_data_froh_l100, aes(x = factor(M), y = FROHP2, fill= factor(M))) +
  geom_boxplot() +
  labs(
    x = "\nMigrants per decade",
    y = "FROH (Crater)\n",
     title = "FROH for Crater at year 2121 for K = 100"
  ) +
    theme_classic()+
    scale_fill_manual(values = c("#FFEDA9","#F4BD19","#ED8000","#692613")) +
    theme(axis.title = element_text(size = 16),
          legend.position = "none",
          plot.title = element_text(size = 16),
          axis.text = element_text(size=16)))
ggsave("popgen_results/froh_2121_K100_GSE.png", dpi = 600,  width = 7.5, height = 6.5)

anova_result <- aov(FROHP2 ~ factor(M), data = plot_data_froh_l100)
summary(anova_result)
TukeyHSD(anova_result)


plot_data_froh_l50 <- all_data %>%
  filter(K == 50,
         Year == 55702)
unique(plot_data_froh_l50$M)

(plot <- ggplot(plot_data_froh_l50, aes(x = factor(M), y = FROHP2, fill = factor(M))) +
  geom_boxplot() +
  labs(
    x = "\nMigrants per decade",
    y = "FROH (Crater)\n",
    title = "FROH for Crater at year 2121 for K = 50"
  ) +
  theme_classic()+
    scale_fill_manual(values = c("#FFEDA9","#F4BD19","#ED8000","#692613")) +
      theme(axis.title = element_text(size = 16),
        legend.position = "none",
        plot.title = element_text(size = 16),
        axis.text = element_text(size=16)))
ggsave("popgen_results/froh_2121_K50_GSE.png", dpi = 600, width = 7.5, height = 6.5)

years_to_use <- c(55603, 55612, 55622, 55632, 55642,
                  55652, 55662, 55672, 55682, 55692, 55702)

timeline_data <- all_data_zoo %>%
  filter(Year %in% years_to_use)
year_labels <- c(
  "55603" = "2022",
  "55612" = "2031",
  "55622" = "2041",
  "55632" = "2051",
  "55642" = "2061",
  "55652" = "2071",
  "55662" = "2081",
  "55672" = "2091",
  "55682" = "2101",
  "55692" = "2111",
  "55702" = "2120"
)
timeline_data_gse <- all_data %>%
  filter(Year %in% years_to_use) %>%
  mutate(
    Year_factor = factor(
      as.character(Year),
      levels = as.character(years_to_use),
      labels = year_labels
    )
  )
levels(timeline_data_gse$Year_factor)

ggplot(timeline_data_gse,
       aes(x = factor(Year_factor), y = FROHP2, fill = factor(M))) +
  geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
  scale_fill_manual(values = c("#FFEDA9","#F4BD19","#ED8000","#692613")) +
  facet_grid(M ~ K, labeller = label_both, scales = "free_y") +
  theme_bw(base_size = 14) +
  labs(
    x = "\nYear",
    y = "FROH\n",
    fill = "Migrants per decade (M)"
  ) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.spacing = unit(1, "lines")
  )
ggsave("popgen_results/froh_gse_evolution.png", dpi = 600, width = 9.5, height = 6.5)

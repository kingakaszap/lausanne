library(tidyverse)

# base folder containing the 8 subfolders
base_path_zoo <- "./popgen_results/Forward_Sim_migrants_from_ZOO/"

# list all csv files inside all subfolders
all_files_zoo <- list.files(base_path_zoo, 
                        pattern = "\\.csv$", 
                        recursive = TRUE, 
                        full.names = TRUE)

# read all files
all_data_zoo <- map_df(all_files_zoo, function(f) {
  
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
head(all_data_zoo)
str(all_data_zoo)
# View(all_data)
str(all_data_zoo)
unique(all_data_zoo$M)
all_data_zoo<-data.frame((all_data_zoo))
all_data_zoo <- all_data_zoo %>%
  mutate(Year = map_int(Year, as.integer))

plot_data_froh_100_zoo <- all_data_zoo %>%
  filter(K == 100,
         Year == 55702)
unique(plot_data_froh_100_zoo$M)

(plot <- ggplot(plot_data_froh_100_zoo, aes(x = factor(M), y = FROHP2, fill = factor(M))) +
    geom_boxplot()+
  labs(
    x = "\nMigrants per decade",
    y = "FROH (Crater)\n",
    title = "FROH for Crater at year 2121 for K = 100"
  ) +
    theme_classic()+
    scale_y_continuous(limits = c(0.20, 0.65)) + 
    
    scale_fill_manual(values = c("#FFEDA9","#F4BD19","#ED8000","#692613")) +
    theme(axis.title = element_text(size = 24),
          legend.position = "none",
          plot.title = element_text(size = 16),
          axis.text = element_text(size=24)))
ggsave("popgen_results/froh_2121_K100_ZOO.png", dpi = 600,  width = 7.5, height = 6.5)
anova_result_zoo <- aov(FROHP2 ~ factor(M), data = plot_data_froh_100_zoo)
summary(anova_result_zoo)
TukeyHSD(anova_result_zoo)


plot_data_froh_50_zoo <- all_data_zoo %>%
  filter(K == 50,
         Year == 55702)
unique(plot_data_froh_50_zoo$M)

(plot <- ggplot(plot_data_froh_50_zoo, aes(x = factor(M), y = FROHP2, fill = factor(M))) +
  geom_boxplot()+
  labs(
    x = "\nMigrants per decade",
    y = "FROH (Crater)\n",
    title = "FROH for Crater at year 2121 for K = 50"
  ) +
    scale_y_continuous(limits = c(0.40, 0.65)) + 
    theme_classic()+
    scale_fill_manual(values = c("#FFEDA9","#F4BD19","#ED8000","#692613")) +
    theme(axis.title = element_text(size = 16),
          legend.position = "none",
          plot.title = element_text(size = 16),
          axis.text = element_text(size=16)))
ggsave("popgen_results/froh_2121_K50_ZOO.png", dpi = 600,  width = 7.5, height = 6.5)

# timeline
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
timeline_data <- all_data_zoo %>%
  filter(Year %in% years_to_use) %>%
  mutate(
    Year_factor = factor(
      as.character(Year),
      levels = as.character(years_to_use),
      labels = year_labels
    )
  )

ggplot(timeline_data,
       aes(x = factor(Year_factor), y = FROHP3, fill = factor(M))) +
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
ggsave("popgen_results/froh_zoo_evolution_FOR_ZOO.png", dpi = 600, width = 9.5, height = 6.5)


# shared plot
plot_data_froh_combined_zoo <- all_data_zoo %>%
  filter(
    K %in% c(50, 100),
    Year == 55702
  )
ggplot(plot_data_froh_combined_zoo,
       aes(x = factor(K), y = FROHP2, fill = factor(M))) +
  
  geom_boxplot() +
  
  scale_fill_manual(values = c("#FFEDA9","#F4BD19","#ED8000","#692613"),
                    name = "Migrants per decade (M)") +
  
  scale_y_continuous(limits = c(0.22, 0.70)) +
  
  labs(
    x = "\nCarrying capacity (K)",
    y = "FROH (Crater)\n",
    title = "Zoo-Crater migration"
  ) +
  
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.position = "none"
  )

ggsave("popgen_results/froh_2121_ALL_ZOO.png", dpi = 600,  width = 7, height = 5)



plot_data_relchange_ZOO <- all_data_zoo %>%
  group_by(r) %>%
  mutate(
    froh_ref = FROHP2[Year == 55604][1],   # value at Year 55604 for this r
    RelFROH = (FROHP2 / froh_ref) * 100 - 100
  ) %>%
  ungroup()

RELATIVE_plot_data_froh_combined_ZOO <- plot_data_relchange_ZOO %>%
  filter(
    K %in% c(50, 100),
    Year == 55702
  )
ggplot(RELATIVE_plot_data_froh_combined_ZOO,
       aes(x = factor(K), y = RelFROH, fill = factor(M))) +
  
  geom_boxplot() +
  
  scale_fill_manual(values = c("#FFEDA9","#F4BD19","#ED8000","#692613"),
                    name = "Migrants per decade (M)") +
  
  # scale_y_continuous(limits = c(0.22, 0.70)) +
  
  labs(
    x = "\nCarrying capacity (K)",
    y = "Change in FROH compared to 2021 (%)\n",
    title = "Zoo-Crater migration"
  ) +
  # geom_line() +
  scale_y_continuous(limits = c(-25, 60)) +
  
  geom_hline(yintercept = 0, linetype = "dashed", lwd = 0.8, colour = "#5E5E5E") +
  theme_classic() +
  theme(
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    plot.title = element_text(size = 16),
    axis.text = element_text(size = 16),
    legend.position = "none"
  )
ggsave("popgen_results/froh_2121_ALL_zoo_relative.png", dpi = 600,  width = 7, height = 5)





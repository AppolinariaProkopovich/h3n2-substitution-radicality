library(readr)
library(dplyr)

df <- read.csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad.csv")
glimpse(df)

# Сайты для исключения
concordant_sites <- c(17, 176, 202, 210, 238, 241, 364, 380, 530)
r_sites <- c(171, 183, 217, 236, 242, 245, 377, 391, 440, 441, 454, 455, 466)
sites_to_exclude <- unique(c(concordant_sites, r_sites))

# Сначала аннотируем ВСЕ данные
df <- df %>%
  mutate(
    gwas_type = case_when(
      site %in% concordant_sites ~ "concordant",
      site %in% r_sites ~ "r_site",
      TRUE ~ "neutral"
    )
  )

# Статистика
cat("Исходные данные:", nrow(df), "строк,", length(unique(df$site)), "сайтов\n")
cat("Concordant строк:", sum(df$gwas_type == "concordant"), "\n")
cat("R-site строк:", sum(df$gwas_type == "r_site"), "\n")
cat("Neutral строк:", sum(df$gwas_type == "neutral"), "\n")

# Датасет БЕЗ артефактных сайтов
df_clean <- df %>% filter(gwas_type == "neutral")
cat("После фильтрации:", nrow(df_clean), "строк,", length(unique(df_clean$site)), "сайтов\n")

write.csv(df_clean, "/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_filtered_concordant_r_site.csv", row.names = FALSE)


library(readr)
library(dplyr)

df <- read.csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_site_Sen_Entr.csv")
glimpse(df)

# Стареющие сайты (senescence, Popova et al. 2019)
sites_to_keep <- c(205, 174, 19, 189, 228, 176, 158, 110, 144, 328, 363, 242)

#  фильтруем
filtered_df <- df %>% filter(site %in% sites_to_keep)
glimpse(filtered_df)


# разметка по классам сайтов, ассоциированных с лабораторным культивированием
# (анализ EpiStat): concordant - замена совпадает с переходом в лабораторные
# условия, r_site - >60% замен на лабораторных ветвях
filtered_df <- filtered_df %>%
  mutate(
    location = case_when(
      site %in% c(17, 176, 202, 210, 238, 241, 364, 380, 530) ~ "concordant",
      site %in% c(171, 183, 217, 236, 242, 245, 377, 391, 440, 441, 454, 455, 466) ~ "r_site",
      TRUE ~ "neutral"
    )
  )
glimpse(filtered_df)

#проверяю что нет unknow
 filtered_df%>% 
  filter(location == "Unknown") %>%
  glimpse()


write.csv(filtered_df, "/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_site_Sen.csv",  row.names=FALSE)
 
model<- lm(miyata_rad ~ distance, data=df)
summary(model)

library(dplyr)

# Для D1 стареющих
d1_sen <- read.csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_site_Sen.csv")
# Для D2 стареющих  
d2_sen <- read.csv("/home/ria/HSE/lab/new_three/data/old_data/old_new_h3_distances_rad_site_sen.csv")

# Для каждого: сводка по distance
summary_table <- bind_rows(
  d1_sen %>% summarise(
    dataset = "D1 (с лаб.)",
    N = n(),
    median_dist = median(distance),
    max_dist = max(distance),
    q75_dist = quantile(distance, 0.75),
    q95_dist = quantile(distance, 0.95)
  ),
  d2_sen %>% summarise(
    dataset = "D2 (клин.)",
    N = n(),
    median_dist = median(distance),
    max_dist = max(distance),
    q75_dist = quantile(distance, 0.75),
    q95_dist = quantile(distance, 0.95)
  )
)

print(summary_table)

head(d1_sen$distance, 10)
head(d2_sen$distance, 10)
# D1 после исключения сайтов, ассоциированных с лабораторным культивированием
d1_filtered <- read.csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_filtered_concordant_r_site_sen.csv")

bind_rows(
  d1_sen %>% summarise(dataset="D1 стареющие", N=n(),
    med=median(distance), q95=quantile(distance,0.95), max=max(distance)),
  d1_filtered %>% summarise(dataset="D1 без ассоциированных сайтов, стареющие", N=n(),
    med=median(distance), q95=quantile(distance,0.95), max=max(distance))
)



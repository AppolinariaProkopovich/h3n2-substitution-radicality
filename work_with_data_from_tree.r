library(readr)
library(dplyr)

df_new <- read_csv("/home/ria/HSE/lab/new_three/data/my_output/h1pdm_test_dif_distsnce2.сsv")
df_new <- df_new %>% select(-der_node)
glimpse(df_new)



df_old <- read_csv("/home/ria/HSE/lab/model/knownpass/h1pdm_distances.csv")
glimpse(df_old)



#сравнить таблицы по колонкам

df_new <- df_new %>% select(-der_node)

df_new <- df_new %>% mutate(across(ends_with("codon"), toupper))
df_old <- df_old %>% mutate(across(ends_with("codon"), toupper))
glimpse(df_old)


key_cols <- c("site", "ancestor_node", "grandanc_codon", "grandanc_aa", 
              "anc_codon", "anc_aa", "der_aa")
only_in_new <- anti_join(df_new, df_old, by=key_cols)
glimpse(only_in_new)

only_in_old <- anti_join(df_old, df_new, by=key_cols)
glimpse(only_in_old)

#--------------------------------
# отличие в distance 
key_cols <- c("site", "grandanc_codon", "grandanc_aa", 
              "anc_codon", "anc_aa", "der_aa")

cmp <- inner_join(
  df_new %>% mutate(across(ends_with("codon"), toupper)),
  df_old %>% mutate(across(ends_with("codon"), toupper)),
  by = key_cols,
  suffix = c("_new", "_old"),
  relationship = "many-to-many"
)

# Отсортируем обе distance внутри каждой группы ключей и сравним
cmp <- cmp %>%
  group_by(across(all_of(key_cols))) %>%
  arrange(distance_old, .by_group = TRUE) %>%
  mutate(rank_old = row_number()) %>%
  ungroup()

cmp <- cmp %>%
  group_by(across(all_of(key_cols))) %>%
  arrange(distance_new, .by_group = TRUE) %>%
  mutate(rank_new = row_number()) %>%
  ungroup()

# Оставляем только пары с совпадающим рангом (т.е. сопоставляем по порядку)
cmp_matched <- cmp %>% filter(rank_old == rank_new)

cat("Сопоставлено строк:", nrow(cmp_matched), "\n")

cmp_matched <- cmp_matched %>%
  mutate(
    dist_diff  = distance_new - distance_old,
    dist_ratio = distance_new / distance_old
  )

cat("\n--- Разница (new - old) ---\n")
summary(cmp_matched$dist_diff)

cat("\n--- Отношение (new / old) ---\n")
summary(cmp_matched$dist_ratio)

# Все ли разницы одинаковы по знаку?
cat("\nВсе отрицательные?", all(cmp_matched$dist_diff < 0), "\n")
cat("Все положительные?", all(cmp_matched$dist_diff > 0), "\n")
cat("Есть нулевые (совпадают)?", any(abs(cmp_matched$dist_diff) < 1e-10), "\n")




#----------------------------------------------------
# смотрела есть ли строки с одинаковой distance и node
df <- read.csv("")

df_df <- df %>%
  group_by(ancestor_node, distance) %>%
  filter(n() > 1) %>%
  arrange(ancestor_node, distance)


glimpse(df_df)

df2 <-  read.csv("")

df_df2 <- df2 %>%
  group_by(ancestor_node, distance) %>%
  filter(n() > 1) %>%
  arrange(ancestor_node, distance)
glimpse(df_df2)


#---------------- new ----------
# Выберем сайт с 5+ заменами для анализа
popular_site <- site_counts$site[site_counts$n >= 5][1]
cat("\nАнализируем сайт:", popular_site, "\n")

site_data <- df_old %>% filter(site == popular_site)
site_data %>% 
  select(ancestor_node, grandanc_aa, anc_aa, der_aa, distance) %>% 
  arrange(ancestor_node, distance)

# Для каждого ancestor_node посмотрим расстояния
cat("\nПо ancestor_node:\n")
for(node in unique(site_data$ancestor_node)) {
  cat("\n", node, ":\n", sep = "")
  subset_data <- site_data %>% filter(ancestor_node == node)
  print(subset_data %>% select(grandanc_aa, anc_aa, der_aa, distance))
}

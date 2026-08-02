library(dplyr)
library(ggplot2)
library(purrr)

# фильтруем такие сайты, в которых от одного узла ancestor_node и от одного сайта
# (соответвенно это конкретная родительская аминокислота) произошли не менее 10 дочерних замен

df <- read.csv("/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_granddist_rad.сsv")

df_filtered <- df %>%
  group_by(ancestor_node, site) %>%
  filter(n() >= 10) %>%
  ungroup()

group_sizes <- df_filtered %>%
  count(ancestor_node, site, name = "n_obs") %>%
  mutate(label = paste0(ancestor_node, " | site=", site)) %>%
  arrange(n_obs) %>%                         # сортировка от меньшего к большему
  mutate(label = factor(label, levels = label)) 

nrow(group_sizes)          
min(group_sizes$n_obs)   


ggplot(group_sizes, aes(x = label, y = n_obs)) +
  geom_col() +
  scale_y_continuous(
    breaks = seq(0, max(group_sizes$n_obs), by = 10),
    minor_breaks = seq(0, max(group_sizes$n_obs), by = 5)
  )+
  labs(
    x = "ancestor_node | site (от меньшего к большему)",
    y = "Число наблюдений",
  ) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
glimpse(df)

#write.csv(df_filtered, file = "/home/ria/HSE/lab/new_three/data/my_output/hipervar_n10.сsv", row.names=FALSE)

  # Посмотри, сколько уникальных значений distance в каждой группе
df_filtered %>%
  group_by(ancestor_node, site) %>%
  summarise(
    unique_miyata = n_distinct(miyata_rad),
    unique_gr     = n_distinct(gr_rad_anc),
    .groups = "drop"
  ) %>%
  filter(unique_miyata == 1 | unique_gr == 1)

  # исключаю  этот сайт

  df_for_lm <- df_filtered %>%
  group_by(ancestor_node, site, anc_aa) %>%
  filter(n_distinct(distance) >= 2) %>%
  filter(n_distinct(miyata_rad) >= 2, n_distinct(gr_rad_anc) >= 2) %>%
  ungroup()

glimpse(df_for_lm)


# Регрессии для каждой группы по gr_rad_anc и miyata_rad

run_regression <- function(data, metric) {
  mod <- lm(reformulate("distance", response = metric), data = data)
  s <- summary(mod)
  coef_row <- s$coefficients["distance", ]
  tibble(
    metric    = metric,
    estimate  = coef_row["Estimate"],
    sign      = ifelse(coef_row["Estimate"] > 0, "+", "-"),
    p_value   = coef_row["Pr(>|t|)"]
  )
}

results <- df_for_lm %>%
  group_by(ancestor_node, site, anc_aa) %>%
  group_modify(~ {
    bind_rows(
      run_regression(.x, "miyata_rad"),
      run_regression(.x, "gr_rad_anc")
    )
  }) %>%
  ungroup()


print(results, n=108)

# если results уже посчитаны:
results_bn <- results %>%
  mutate(p_bonferroni = p.adjust(p_value, method = "bonferroni"),
         p_fdr        = p.adjust(p_value, method = "fdr"))

# Бонферрони — строгая поправка
results_bn %>% filter(p_bonferroni < 0.05)
print(results_bn, n=108)

# FDR (Benjamini-Hochberg) — мягче, контролирует долю ложных среди значимых
results_bn %>% filter(p_fdr < 0.05)


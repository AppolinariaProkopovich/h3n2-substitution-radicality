library(tidyverse)
library(broom)

# ══════════════════════════════════════════════════════════════════════════
# 1. Загрузка данных
# ══════════════════════════════════════════════════════════════════════════

df <- read.csv("/home/ria/HSE/lab/new_three/data/old_data/old_new_h3_distances_rad.csv")
glimpse(df)

# Колонки в датасете:
#   gr_rad_anc        = Grantham(anc→der)
#   gr_rad_anc_grand  = Grantham(grandanc→der)
#   miyata_rad        = Miyata(anc→der)
#   miyata_rad_grand  = Miyata(grandanc→der)
#   distance          = филогенетическое расстояние anc→der
#   distance_grand    = филогенетическое расстояние grandanc→der

# ══════════════════════════════════════════════════════════════════════════
# 2. Разметка сайтов (бинарные колонки — сайт может быть в нескольких)
# ══════════════════════════════════════════════════════════════════════════

h3_surface <- c(25,37,38,41,43,47,48,49,54,61,62,64,66,71,73,79,94,97,98,99,
                107,108,110,112,120,140,142,144,147,148,149,151,153,156,158,
                159,160,161,173,174,175,176,178,179,181,183,187,188,189,204,
                205,206,208,209,214,215,224,228,230,238,240,241,255,256,277,
                278,279,280,285,287,289,292,294,301,305,307,326,328,329,340,
                341,352,353,356,361,363,364,372,374,375,376,377,383,384,387,
                391,394,398,402,403,405,406,414,416,418,427,466,480,488,491,
                492,499,500,501,503,505,506,509,510,513,517,518)

h3_internal <- c(27,29,31,32,33,35,42,44,45,50,52,53,58,59,60,67,68,72,75,77,
                 80,82,83,84,85,86,87,88,89,92,95,100,102,103,104,105,106,113,
                 114,115,118,123,124,125,126,127,128,129,131,132,133,134,136,
                 141,143,146,155,163,164,166,167,168,169,170,177,180,182,186,
                 192,193,194,195,196,197,198,199,200,201,202,207,211,216,218,
                 219,220,221,222,229,231,236,244,245,246,247,248,251,253,257,
                 259,260,261,263,265,266,267,268,269,270,272,273,274,281,282,
                 283,284,286,288,291,297,298,299,302,303,304,310,311,318,319,
                 321,322,325,330,332,333,335,336,338,349,350,351,355,358,359,
                 362,367,368,369,373,385,386,389,393,396,411,420,423,425,426,
                 428,429,430,432,434,435,436,437,438,440,441,444,445,446,448,
                 449,452,453,454,455,456,457,458,459,460,463,464,467,469,471,
                 474,476,477,481,483,485,486,487,489,493,494,497,502,508,511,
                 512,515)

entr_sites  <- c(69, 291, 206, 111)
sen_sites   <- c(205, 174, 19, 189, 228, 176, 158, 110, 144, 328, 363, 242)
clust_sites <- c(144, 66, 204, 277, 22, 305, 388, 13, 184, 209, 158, 509,
                 3, 109, 281, 422, 98, 521, 45, 313, 128)

df <- df %>%
  mutate(
    is_surface  = as.integer(site %in% h3_surface),
    is_internal = as.integer(site %in% h3_internal),
    is_entr     = as.integer(site %in% entr_sites),
    is_sen      = as.integer(site %in% sen_sites),
    is_clust    = as.integer(site %in% clust_sites)
  )

# ══════════════════════════════════════════════════════════════════════════
# 3. Функция регрессии + логирование пропусков
# ══════════════════════════════════════════════════════════════════════════

# Глобальный список для пропущенных групп
skipped_groups <- list()

safe_lm <- function(data, y_col, x_col, site_val, anc_aa_val, metric, predictor) {
  n_obs       <- nrow(data)
  n_distinct_x <- n_distinct(data[[x_col]])
  
  # Причина пропуска
  reason <- NA_character_
  if (n_obs < 3) {
    reason <- paste0("мало наблюдений (n=", n_obs, ")")
  } else if (n_distinct_x < 2) {
    reason <- paste0("одинаковый предиктор (", x_col, "=", round(data[[x_col]][1], 4), ", n=", n_obs, ")")
  }
  
  if (!is.na(reason)) {
    # Записываем в лог пропусков
    skipped_groups[[length(skipped_groups) + 1]] <<- tibble(
      site = site_val, anc_aa = anc_aa_val,
      metric = metric, predictor = predictor,
      n_obs = n_obs, reason = reason
    )
    return(tibble(
      estimate = NA_real_, std.error = NA_real_,
      statistic = NA_real_, p.value = NA_real_,
      r.squared = NA_real_, n_obs = n_obs
    ))
  }
  
  mod <- lm(reformulate(x_col, y_col), data = data)
  tidy(mod) %>%
    filter(term == x_col) %>%
    select(estimate, std.error, statistic, p.value) %>%
    mutate(r.squared = summary(mod)$r.squared, n_obs = n_obs)
}

# ══════════════════════════════════════════════════════════════════════════
# 4. Перебор всех комбинаций
# ══════════════════════════════════════════════════════════════════════════

#combos <- tribble(
 # ~y_col,              ~x_col,          ~metric,    ~predictor,
  #"miyata_rad",        "distance",       "miyata",   "distance",
  #"gr_rad_anc",        "distance",       "grantham", "distance",
  #"miyata_rad_grand",  "distance_grand", "miyata",   "granddistance",
  #"gr_rad_anc_grand",  "distance_grand", "grantham", "granddistance",
#)

combos <- tribble(
  ~y_col,       ~x_col,     ~metric,    ~predictor,
  "miyata_rad",  "distance", "miyata",   "distance",
  "gr_rad_anc",  "distance", "grantham", "distance",
)

# Группируем по site × anc_aa
nested <- df %>%
  group_by(site, anc_aa) %>%
  group_nest()

# Все комбинации подгрупп × метрик
all_tasks <- nested %>% cross_join(combos)

# Запускаем регрессии
results <- all_tasks %>%
  rowwise() %>%
  mutate(
    res = list(safe_lm(data, y_col, x_col, site, anc_aa, metric, predictor))
  ) %>%
  unnest(res) %>%
  select(site, anc_aa, metric, predictor,
         estimate, std.error, statistic, p.value,
         r.squared, n_obs) %>%
  ungroup()

# Добавляем метки сайтов
site_labels <- df %>%
  distinct(site, is_surface, is_internal, is_entr, is_sen, is_clust)

results <- results %>%
  left_join(site_labels, by = "site") %>%
  arrange(site, anc_aa, metric, predictor)

# Разделяем на результаты и пропуски
results_clean  <- results %>% filter(!is.na(estimate))
skipped_df     <- bind_rows(skipped_groups)

# ══════════════════════════════════════════════════════════════════════════
# 5. Сохраняем
# ══════════════════════════════════════════════════════════════════════════

#write_tsv(results_clean, "/home/ria/HSE/lab/new_three/data/old_data/old_ new_data_site_anc_aa/old_new_data_grouped_regressions_results.tsv")
#write_tsv(skipped_df, "/home/ria/HSE/lab/new_three/data/old_data/old_ new_data_site_anc_aa/old_new_data_grouped_regressions_skipped.tsv")

cat("=== Результаты ===\n")
cat("Подгрупп site x anc_aa:", nrow(nested), "\n")
cat("Задач (подгрупп x 4 комбинации):", nrow(all_tasks), "\n")
cat("Успешных регрессий:", nrow(results_clean), "\n")
cat("Пропущено:", nrow(skipped_df), "\n\n")

# ══════════════════════════════════════════════════════════════════════════
# 6. Лог пропущенных групп
# ══════════════════════════════════════════════════════════════════════════

cat("=== Пропущенные группы (первые 20) ===\n")
skipped_df %>%
  print(n = 20)

cat("\n=== Причины пропусков ===\n")
skipped_df %>%
  mutate(reason_type = ifelse(str_detect(reason, "мало"), "мало наблюдений", "одинаковый предиктор")) %>%
  count(reason_type) %>%
  print()

# ══════════════════════════════════════════════════════════════════════════
# 7. Сводки
# ══════════════════════════════════════════════════════════════════════════

cat("\n=== Сводка по метрика x предиктор ===\n")
results_clean %>%
  group_by(metric, predictor) %>%
  summarise(
    n_groups     = n(),
    n_signif     = sum(p.value < 0.05, na.rm = TRUE),
    pct_negative = round(mean(estimate < 0, na.rm = TRUE) * 100, 1),
    median_r2    = round(median(r.squared, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print()

cat("\n=== Сводка по биологическим группам (Sen / Entr / Clust) ===\n")
results_clean %>%
  pivot_longer(cols = c(is_entr, is_sen, is_clust),
               names_to = "group", values_to = "flag") %>%
  filter(flag == 1) %>%
  mutate(group = str_remove(group, "is_")) %>%
  group_by(group, metric, predictor) %>%
  summarise(
    n_groups   = n(),
    n_signif   = sum(p.value < 0.05),
    median_est = round(median(estimate, na.rm = TRUE), 4),
    median_r2  = round(median(r.squared, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print(n = 40)

cat("\n=== Сводка по расположению (surface / internal) ===\n")
results_clean %>%
  pivot_longer(cols = c(is_surface, is_internal),
               names_to = "location", values_to = "flag") %>%
  filter(flag == 1) %>%
  mutate(location = str_remove(location, "is_")) %>%
  group_by(location, metric, predictor) %>%
  summarise(
    n_groups   = n(),
    n_signif   = sum(p.value < 0.05, na.rm = TRUE),
    median_est = round(median(estimate, na.rm = TRUE), 4),
    median_r2  = round(median(r.squared, na.rm = TRUE), 3),
    .groups = "drop"
  ) %>%
  print(n = 20)

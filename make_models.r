library(dplyr)
library(lmerTest)
library(readr)
library(ggplot2)

packageVersion("lmerTest")
packageVersion("lme4")
packageVersion("rlang")


"ancestor_node - Узел в филогенетическом дереве, где произошла замена.
grandanc_codon - Кодон-предок (до предыдущей замены).
grandanc_aa -
anc_codon - Кодон-родитель (результат предыдущей замены)
anc_aa - 
der_codon - кодон (результат текущей замены).
der_aa
distance - 
Sen - senescence - стареющие
Entr - entrenchment - акклиматизирующиеся
 "

#группируем
group_counts <- datasets_rad %>%
  group_by(site) %>%
  tally() 

#фильтруем
valid_groups <- group_counts %>%
  filter(n >= 3) %>%
  pull(site)

#создаем отфильтрованый dataset
datasets_rad_site <- datasets_rad %>%
  filter(anc_aa %in% valid_groups)
glimpse(datasets_rad_site)

# min группа
group_min <- datasets_rad_site %>%
  group_by(site) %>%
  tally()
 min_n <- min(group_min$n)
print(min_n)

#n групп
group_n <- datasets_rad_site %>% group_by(site)
glimpse(group_n)
group_n

#для вложенных
# фильтруем  сайты в которых не менее 3 разных anc_aa
datasets_rad_site_anc_aa <- datasets_rad_site %>%
  group_by(site) %>%
  filter(n_distinct(anc_aa) == 3) %>%
  ungroup()

glimpse(datasets_rad_site_anc_aa)

  datasets_rad_site_anc_aa %>%
  group_by(site) %>%
  summarise(
    n_anc_aa = n_distinct(anc_aa),  # подгрупп в сайте
    n_obs = n()                      # наблюдений в сайте
  ) %>%
  summary()
 


#write.csv(datasets_rad_site, '/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist_site_anc_aa.сsv', row.names=FALSE)


datasets_rad <-  read_csv('/home/ria/HSE/lab/new_three/data/old_data/clust/old_data_h3_distances_rad_clust_509.csv')
glimpse(datasets_rad)


#MОДЕЛИ 
model<- lm(gr_rad_anc ~ distance, data=datasets_rad)
summary(model)

model<- lm(miyata_rad ~ distance, data=datasets_rad)
summary(model)

model <- lmer(miyata_rad ~ distance  + (1|site), data=datasets_rad_site)
summary(model)

model <- lmer(gr_rad_anc ~ distance  + (1|site), data=datasets_rad_site_anc_aa)
summary(model)

model <- lmer(gr_rad_anc  ~ distance + (distance|anc_aa), data=datasets_rad_serf_inter_site ) 
summary(model)

model <- lmer(gr_rad_anc  ~ distance  + (distance||grandanc_aa), data=datasets_rad_serf_inter_site ) 
summary(model)

model <- lmer(gr_rad_anc  ~ distance + (distance|der_aa) , data=datasets_rad_serf_inter_site ) 
summary(model)

model <- lmer(miyata_rad ~ distance + (distance||location), data=datasets_rad_serf_inter_site)
summary(model)

model_site <- lmer(miyata_rad ~ distance  + (1|site), data=datasets_rad_site)
summary(model_site)

model_site_anc_aa <- lmer(miyata_rad ~ distance + (1|site/anc_aa), data=datasets_rad_site)
summary(model_site_anc_aa)

AIC(model_site, model_site_anc_aa)

# Если в датасете есть колонка с типом белка:
p <- ggplot(datasets_rad, aes(x = distance, y = miyata_rad, color = factor(site))) +
  geom_point(alpha = 0.3) +
  theme_minimal() +
  theme(legend.position = "none")

print(p)

"разные сайты покрывают разные диапазоны distance. 
Некоторые цвета присутствуют только на коротких дистанциях (0–5), другие — только на длинных (10–18). 
Именно это создаёт ложный тренд в модели без случайных эффектов: если сайты с высокой радикальностью 
случайно чаще встречаются при коротких distance, а сайты с низкой радикальностью — при длинных, 
линейная модель «видит» снижение радикальности со временем, хотя на самом деле это просто разные сайты."

p2 <- ggplot(datasets_rad, aes(x = distance)) +
  geom_histogram(bins = 30, fill = "steelblue") +
  facet_wrap(~site, scales = "free_y") +
  theme_minimal() +
  theme(strip.text = element_text(size = 6))

print(p2)

" Разные сайты живут в разных диапазонах distance и имеют разные уровни радикальности. 
Без учёта сайта модель путает «какой это сайт» с «сколько прошло времени». 
Именно поэтому наклон исчезает при добавлении (1|site)"


datasets_rad_site$distance_sc <- scale(datasets_rad_site$distance, center = TRUE, scale = TRUE)

model <- lmer(gr_rad_anc ~ distance_sc + (distance_sc | site),
              data = datasets_rad_site)
summary(model)



model_bobyqa <- lmer(
  gr_rad_anc ~ distance_sc + (distance_sc | site),
  data = datasets_rad_site,
  control = lmerControl(optimizer = "bobyqa",
                        optCtrl = list(maxfun = 2e5))
)


isSingular(model_bobyqa, tol = 1e-5)



VarCorr(model_bobyqa)



# отдельно стареющие и акклиматизирующиеся
ds_s <- datasets_rad_sen_entr %>% filter(selection_type == "Sen")
ds_e <- df%>% filter(selection_type == "Entr")

model<- lm(gr_rad_anc ~ distance , data=ds_s)
summary(model)
model<- lmer(gr_rad_anc ~ distance + (1|site), data=ds_e)
summary(model)



# отдельно serface и internal
datasets_rad <- read.csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist_serf_inter.csv')
glimpse(datasets_rad)

ds_s <- datasets_rad %>% filter(location == "serf")
ds_i <- datasets_rad %>% filter(location == "inter")


model<- lm(gr_rad_anc ~ distance, data=ds_s)
summary(model)
model<- lm(miyata_rad ~ distance, data=ds_s)
summary(model)

model<- lm(gr_rad_anc ~ distance, data=ds_i)
summary(model)
model<- lm(miyata_rad ~ distance, data=ds_i)
summary(model)



# визуалтзация данных
# Сколько наблюдений в каждой группе
table(datasets_rad$location)

# вариация внутри групп
tapply(datasets_rad$gr_rad_anc, 
       datasets_rad$location, 
       var)


ggplot(datasets_rad_sen_entr, aes(x=distance, y=miyata_rad, color=selection_type)) +
  geom_point(alpha=0.5) +
  geom_smooth(method="lm", se=FALSE) +
  facet_wrap(~selection_type) +
  theme_bw()

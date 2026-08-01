# Загрузить данные D1 (все замены с сайтами и радикальностью)
d <- read.csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad.csv")

# Список 24 GWAS-сайтов
gwas_sites <- c(17, 176, 202, 210, 238, 241, 364, 380, 530,  # concordant
                171, 183, 217, 236, 242, 245, 377, 391, 
                440, 441, 454, 455, 466, 468, 558)           # R-sites

# Пометить замены
d$is_gwas <- d$site %in% gwas_sites

# Сравнение радикальности: GWAS-сайты vs остальные
# (с учётом родительской аминокислоты)
wilcox.test(miyata_rad ~ is_gwas, data = d)
wilcox.test(gr_rad_anc ~ is_gwas, data = d)

# Средние
aggregate(cbind(miyata_rad, gr_rad_anc) ~ is_gwas, data = d, 
          FUN = function(x) c(mean = mean(x), median = median(x), n = length(x)))

# Более строгое сравнение: для каждой родительской аминокислоты
# (контроль за составом замен)
library(lme4)
library(lmerTest)
m_miyata <- lmer(miyata_rad ~ is_gwas + (1|anc_aa), data = d)
summary(m_miyata)

m_grantham <- lmer(gr_rad_anc ~ is_gwas + (1|anc_aa), data = d)
summary(m_grantham)

# Или ещё строже: с контролем за сайтом И родительской аминокислотой
# (сравниваем замены из GWAS-сайтов с заменами из тех же anc_aa в других сайтах)
m2 <- lm(miyata_rad ~ is_gwas + anc_aa, data = d)
summary(m2)


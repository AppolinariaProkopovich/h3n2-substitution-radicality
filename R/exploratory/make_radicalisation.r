library(dplyr)
library(stringr)
library(grantham)



datasets <- read.csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_granddistance/old_data_h3_granddistances.csv")
glimpse(datasets)


# GRATHAM

aa_code <- data.frame(
  one_letter = c("A","R","N","D","C","Q","E","G","H",
                 "I","L","K","M","F","P","S","T","W","Y","V"),
  three_letter = c("Ala","Arg","Asn","Asp","Cys","Gln","Glu","Gly","His",
                   "Ile","Leu","Lys","Met","Phe","Pro","Ser","Thr","Trp","Tyr","Val")
)

allowed <- aa_code$one_letter

# Убираем строки с недопустимыми аминокислотами
bad_letters <- setdiff(unique(c(datasets$anc_aa, datasets$der_aa)), allowed)
bad_letters  # посмотрите что выкидываете

datasets <- datasets %>%
  filter(
    anc_aa %in% allowed,
    der_aa %in% allowed,
    grandanc_aa %in% allowed
  )

# Конвертация в трёхбуквенные коды
convert <- function(aa_vector) {
  conversion_vec <- setNames(aa_code$three_letter, aa_code$one_letter)
  dplyr::recode(aa_vector, !!!conversion_vec)
}

datasets$grandanc_aa_tr <- convert(datasets$grandanc_aa)
datasets$anc_aa_tr      <- convert(datasets$anc_aa)
datasets$der_aa_tr      <- convert(datasets$der_aa)

# Радикальность по Grantham (замена anc -> der)
gr_rad <- grantham_distance(datasets$anc_aa_tr, datasets$der_aa_tr) 
datasets$gr_rad_anc <- gr_rad[[3]] 

glimpse(datasets)

# MIYATA 
miyata <- readRDS("/home/ria/HSE/lab/my_model/miyata.rds")

miyata_distance <- function(aa1, aa2) {
  miyata[aa1, aa2]
}

miyata_rad <- mapply(miyata_distance, datasets$grandanc_aa, datasets$der_aa)

datasets_rad <- datasets %>%
  mutate(
    miyata_rad_grand
  )
glimpse(datasets_rad)


#BLOSUM62
# Чем ниже (отрицательнее) score, тем радикальнее замена.

blosum62 <- matrix(c(
# A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V
  4,-1,-2,-2, 0,-1,-1, 0,-2,-1,-1,-1,-1,-2,-1, 1, 0,-3,-2, 0,  # A
 -1, 5, 0,-2,-3, 1, 0,-2, 0,-3,-2, 2,-1,-3,-2,-1,-1,-3,-2,-3,  # R
 -2, 0, 6, 1,-3, 0, 0, 0, 1,-3,-3, 0,-2,-3,-2, 1, 0,-4,-2,-3,  # N
 -2,-2, 1, 6,-3, 0, 2,-1,-1,-3,-4,-1,-3,-3,-1, 0,-1,-4,-3,-3,  # D
  0,-3,-3,-3, 9,-3,-4,-3,-3,-1,-1,-3,-1,-2,-3,-1,-1,-2,-2,-1,  # C
 -1, 1, 0, 0,-3, 5, 2,-2, 0,-3,-2, 1, 0,-3,-1, 0,-1,-2,-1,-2,  # Q
 -1, 0, 0, 2,-4, 2, 5,-2, 0,-3,-3, 1,-2,-3,-1, 0,-1,-3,-2,-2,  # E
  0,-2, 0,-1,-3,-2,-2, 6,-2,-4,-4,-2,-3,-3,-2, 0,-2,-2,-3,-3,  # G
 -2, 0, 1,-1,-3, 0, 0,-2, 8,-3,-3,-1,-2,-1,-2,-1,-2,-2, 2,-3,  # H
 -1,-3,-3,-3,-1,-3,-3,-4,-3, 4, 2,-3, 1, 0,-3,-2,-1,-3,-1, 3,  # I
 -1,-2,-3,-4,-1,-2,-3,-4,-3, 2, 4,-2, 2, 0,-3,-2,-1,-2,-1, 1,  # L
 -1, 2, 0,-1,-3, 1, 1,-2,-1,-3,-2, 5,-1,-3,-1, 0,-1,-3,-2,-2,  # K
 -1,-1,-2,-3,-1, 0,-2,-3,-2, 1, 2,-1, 5, 0,-2,-1,-1,-1,-1, 1,  # M
 -2,-3,-3,-3,-2,-3,-3,-3,-1, 0, 0,-3, 0, 6,-4,-2,-2, 1, 3,-1,  # F
 -1,-2,-2,-1,-3,-1,-1,-2,-2,-3,-3,-1,-2,-4, 7,-1,-1,-4,-3,-2,  # P
  1,-1, 1, 0,-1, 0, 0, 0,-1,-2,-2, 0,-1,-2,-1, 4, 1,-3,-2,-2,  # S
  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 1, 5,-2,-2, 0,  # T
 -3,-3,-4,-4,-2,-2,-3,-2,-2,-3,-2,-3,-1, 1,-4,-3,-2,11, 2,-3,  # W
 -2,-2,-2,-3,-2,-1,-2,-3, 2,-1,-1,-2,-1, 3,-3,-2,-2, 2, 7,-1,  # Y
  0,-3,-3,-3,-1,-2,-2,-3,-3, 3, 1,-2, 1,-1,-2,-2, 0,-3,-1, 4   # V
), nrow = 20, byrow = TRUE)

#подписываем строки и стобцы матрицы
aa_order <- c("A","R","N","D","C","Q","E","G","H",
              "I","L","K","M","F","P","S","T","W","Y","V")
rownames(blosum62) <- aa_order
colnames(blosum62) <- aa_order

#считаем 

datasets_rad$blosum62_score <- mapply(
  function(a1, a2) blosum62[a1, a2],
  datasets_rad$anc_aa, datasets_rad$der_aa
)

glimpse(datasets_rad)

# Z-SCALES 

# Z1 ~ lipophilicity, Z2 ~ steric/size, Z3 ~ polarity,
# Z4 ~ electronic, Z5 ~ electronegativity/heat of formation
# Радикальность = евклидово расстояние в z-пространстве

zscales <- data.frame(
  aa = c("A","R","N","D","C","Q","E","G","H",
         "I","L","K","M","F","P","S","T","W","Y","V"),
  z1 = c( 0.24, 3.52, 3.05, 3.98,-0.26, 1.75, 3.11, 2.05, 2.47, #липофильность
         -3.89,-4.28, 2.29,-2.85,-4.22, 0.72, 2.39, 0.75,-4.36,-2.54,-2.59), 
  z2 = c(-2.32, 2.50,-3.65,-3.80, 0.41,-1.30,-0.11,-4.06, 1.95, # размер
         -1.73,-1.30, 0.89,-1.66, 1.94,-1.07,-1.07, 0.00, 3.94, 2.44,-2.64),
  z3 = c( 0.60,-3.50, 1.09, 1.70, 3.50, 1.30, 0.42,-0.82,-0.27, # полярность
          1.71,-1.49, 0.32, 1.18, 0.56, 1.78,-1.39,-2.18, 0.59, 0.43,-1.54),
  z4 = c(-0.14, 1.45, 0.42,-0.54,-1.34,-0.93,-1.75, 0.22, 0.28, # электронные свойства
          1.04, 0.05,-0.89, 0.18, 1.30,-1.55, 0.67,-1.42, 0.37,-0.63, 1.01),
  z5 = c(-0.30, 0.22,-1.15,-1.46, 0.21, 0.36,-0.07, 1.65, 0.26, # электроотрицательность
         -0.11,-0.11, 0.15, 0.47, 0.78, 0.85,-0.43, 0.08, 0.30, 0.30,-0.25)
)
rownames(zscales) <- zscales$aa

# Евклидово расстояние по 5 z-шкалам
datasets_rad$zscale_dist <- mapply(
  function(a1, a2) {
    z1 <- as.numeric(zscales[a1, 2:6])
    z2 <- as.numeric(zscales[a2, 2:6])
    sqrt(sum((z1 - z2)^2))
  },
  datasets_rad$anc_aa, datasets_rad$der_aa
)

glimpse(datasets_rad)

write_csv(datasets_rad, "/home/ria/HSE/lab/new_three/data/old_data/old_data_granddistance/old_data_h3_granddistances_rad.csv")

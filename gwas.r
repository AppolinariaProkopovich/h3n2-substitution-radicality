library(ape)
install.packages('ape')
library(castor)

tree <- read.tree("/home/ria/HSE/lab/new_three/data/old_data/knownpass/h3.l.r.newick")

# Загрузить категоризацию
cats <- read.delim("/home/ria/HSE/lab/new_three/data/old_data/gwas/final_categorized_strains.tsv")

# Для каждого листа: 1 = CLINICAL, 2 = LABORATORY
tip_states <- ifelse(cats$laboratory_general == "CLINICAL", 1, 2)
# (убедиться, что порядок совпадает с tree$tip.label)

asr <- asr_max_parsimony(tree, tip_states)
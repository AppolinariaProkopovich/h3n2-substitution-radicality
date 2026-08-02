library(readr)

#старые данные - еще с яйцами
datasets_rad <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad.csv")
glimpse(datasets_rad)

datasets_rad_site <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_site.сsv')
glimpse(datasets_rad_site)

datasets_rad_site_anc_aa <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_site_anc_aa.сsv')
glimpse(datasets_rad_site_anc_aa)

# старые данные кластеризующиеся
datasets_rad <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_clust.csv")
glimpse(datasets_rad)

datasets_rad_site <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_clust_site.csv')
glimpse(datasets_rad_site)

datasets_rad_site_anc_aa <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_clust_site_anc_aa.csv')
glimpse(datasets_rad_site_anc_aa)

# старые данные стареющие и акклиматизирующиеся
datasets_rad <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_Sen_Entr.csv")
glimpse(datasets_rad)

datasets_rad_site <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_site_Sen_Entr.csv')
glimpse(datasets_rad_site)

datasets_rad_site_anc_aa <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_site_anc_aa_Sen_Entr.csv')
glimpse(datasets_rad_site_anc_aa)

datasets_rad_site_anc_aa <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_site_anc_aa_Sen.csv')
glimpse(datasets_rad_site_anc_aa)

"/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_site_Sen.csv"


#старые данные - после удаления сайтов concordance и r-sites
datasets_rad <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_filtered_concordant_r_site.csv")
glimpse(datasets_rad)

datasets_rad_site <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_filtered_concordant_r_site_site.csv")
glimpse(datasets_rad_site)

datasets_rad_site_anc_aa <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_filtered_concordant_r_site_site_anc_aa.csv")
glimpse(datasets_rad_site_anc_aa)

# старые данные - после удаления сайтов concordance и r-sites - кластеризующиеся
datasets_rad <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_filtered_concordant_r_site_clust.csv")
glimpse(datasets_rad)

datasets_rad_site <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_filtered_concordant_r_site_site_clust.csv")
glimpse(datasets_rad_site)

datasets_rad_site_anc_aa <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_filtered_concordant_r_site_site_anc_aa_clust.csv")
glimpse(datasets_rad_site_anc_aa)

#старые данные - после удаления сайтов concordance и r-sites - стареющие
datasets_rad <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_filtered_concordant_r_site_sen.csv")
glimpse(datasets_rad)

datasets_rad_site <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_filtered_concordant_r_site_site_sen.csv")
glimpse(datasets_rad_site)

datasets_rad_site_anc_aa <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_data_h3_distances_rad_filtered_concordant_r_site_site_anc_aa_sen.csv")
glimpse(datasets_rad_site_anc_aa)





# old_new - из данных Ники
datasets_rad <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_new_h3_distances_rad.csv")
glimpse(datasets_rad)

datasets_rad_site <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_new_h3_distances_rad_site.csv')
glimpse(datasets_rad_site)

datasets_rad_site_anc_aa <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_new_h3_distances_rad_site_anc_aa.csv')
glimpse(datasets_rad_site_anc_aa)

# old_new - из данных Ники - стареющие
datasets_rad <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_new_h3_distances_rad_sen.csv")
glimpse(datasets_rad)

datasets_rad_site <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_new_h3_distances_rad_site_sen.csv')
glimpse(datasets_rad_site)

datasets_rad_site_anc_aa <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_new_h3_distances_rad_site_anc_aa_sen.csv')
glimpse(datasets_rad_site_anc_aa)


# old_new - из данных Ники - кластеризующиеся

datasets_rad <- read_csv("/home/ria/HSE/lab/new_three/data/old_data/old_new_h3_distances_rad_clust.csv")
glimpse(datasets_rad)

datasets_rad_site <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_new_h3_distances_rad_clust_site.csv')
glimpse(datasets_rad_site)

datasets_rad_site_anc_aa <- read_csv('/home/ria/HSE/lab/new_three/data/old_data/old_new_h3_distances_rad_clust_site_anc_aa.csv')
glimpse(datasets_rad_site_anc_aa)


# новые данные стареющие и акклиматизирующиеся

datasets_rad <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_SEN_entr.сsv')
glimpse(datasets_rad)

datasets_rad <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist_sen.сsv')
glimpse(datasets_rad)

datasets_rad_site_anc_aa <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist_site_anc_aa_sen.сsv')
glimpse(datasets_rad_site_anc_aa)

# новые данные кластеризующиеся

datasets_rad <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist_clust.сsv')
glimpse(datasets_rad)

datasets_rad_site <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist_clust_site.сsv')
glimpse(datasets_rad_site)

datasets_rad_site_anc_aa <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist_clust_site_anc_aa.сsv')
glimpse(datasets_rad_site_anc_aa)


# новые данные
datasets_rad <- read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist.сsv')
glimpse(datasets_rad)

datasets_rad_site <- read_csv("/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist_site.сsv")
glimpse(datasets_rad)


"/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist_site_anc_aa.сsv"


datasets_rad_anc_aa <- read_csv('/home/ria/HSE/lab/new_three/data/my_output/rad/datasets_rad_anc_aa.csv')
glimpse(datasets_rad_anc_aa)

datasets_rad_grandanc_aa <- read_csv('/home/ria/HSE/lab/new_three/data/my_output/rad/datasets_rad_grandanc_aa.csv')
glimpse(datasets_rad_grandanc_aa)

datasets_rad_der_aa <- read_csv('/home/ria/HSE/lab/new_three/data/my_output/rad/datasets_rad_der_aa.csv')
glimpse(datasets_rad_der_aa)

datasets_rad_site_anc_aa <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/rad/datasets_rad_site_anc_aa.csv')
glimpse(datasets_rad_site_anc_aa)

# новые данные поверхностные и внутренние
datasets_rad_serf_inter <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist_serf_inter_site.csv')
glimpse(datasets_rad_serf_inter)


datasets_rad_serf_site <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist_serf_site.csv')
glimpse(datasets_rad_serf_site)


datasets_rad_inter_site <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_dist_inter_site.csv')
glimpse(datasets_rad_inter_site)


datasets_rad_grand_serf_inter_site <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_granddist_rad_serf_inter_site.сsv')
glimpse(datasets_rad_grand_serf_inter_site)

datasets_rad_grand_serf_site <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_granddist_rad_serf_site.сsv')
glimpse(datasets_rad_grand_serf_site)

datasets_rad_grand_inter <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_granddist_rad_iner.сsv')
glimpse(datasets_rad_grand_inter)


datasets_rad_grand_Entr_Sen <-  read_csv('/home/ria/HSE/lab/new_three/data/my_output/treefile_without_outliers_granddist_rad_Entr_Sen.сsv')
glimpse(datasets_rad_grand_Entr_Sen)




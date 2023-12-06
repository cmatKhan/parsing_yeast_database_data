library(tidyverse)
library(readxl)
library(here)

genes_table = read_csv(here("data/genes_table.csv.gz"))

harbison_pval_df_list = list(
  ypd = read_excel(here('data/harbison/pvalbygene_forpaper_abbr.xls'),
                   sheet='YPD', skip=1) %>%
    dplyr::rename(locus_tag=`...1`, gene_name=`...2`, description=`...3`) %>%
    mutate(across(!c(locus_tag,gene_name,description), as.numeric)) %>%
    pivot_longer(!c(locus_tag,gene_name,description), names_to='tf_cond', values_to='pval') %>%
    separate(tf_cond, c('tf','cond'),sep="_"),
  other_conds = read_excel(here('data/harbison/pvalbygene_forpaper_abbr.xls'),
                         sheet = "other conditions", skip=1) %>%
    dplyr::rename(locus_tag=`...1`, gene_name=`...2`, description=`...3`) %>%
    mutate(across(!c(locus_tag,gene_name,description), as.numeric)) %>%
    pivot_longer(!c(locus_tag,gene_name,description), names_to='tf_cond', values_to='pval') %>%
    separate(tf_cond, c('tf','cond'),sep="_")
)

harbison_pval_df = bind_rows(harbison_pval_df_list) %>%
  filter(!is.na(pval))

harbison_binding_ratio_df_list = list(
  ypd = read_excel(here('data/harbison/ratiobygene_forpaper_abbr.xls'),
                   sheet='YPD',skip=1) %>%
    dplyr::rename(locus_tag=`...1`, gene_name=`...2`, description=`...3`) %>%
    mutate(across(!c(locus_tag,gene_name,description), as.numeric)) %>%
    pivot_longer(!c(locus_tag,gene_name,description), names_to='tf_cond', values_to='binding_ratio') %>%
    separate(tf_cond, c('tf','cond'),sep="_"),
  other_conds = read_excel(here('data/harbison/ratiobygene_forpaper_abbr.xls'),
                         sheet='Other Conditions', skip=1) %>%
    dplyr::rename(locus_tag=`...1`, gene_name=`...2`, description=`...3`) %>%
    mutate(across(!c(locus_tag,gene_name,description), as.numeric)) %>%
    pivot_longer(!c(locus_tag,gene_name,description), names_to='tf_cond', values_to='binding_ratio') %>%
    separate(tf_cond, c('tf','cond'),sep="_")
)

harbison_binding_ratio_df = bind_rows(harbison_binding_ratio_df_list) %>%
  filter(!is.na(binding_ratio))

combined_harbison = harbison_binding_ratio_df %>%
  select(-gene_name) %>%
  left_join(select(harbison_pval_df, -gene_name)) %>%
  filter(!is.na(pval)) %>%
  filter(tf != 'A1 (MATA1)')

locus_tags = combined_harbison %>%
  select(locus_tag) %>%
  distinct() %>%
  left_join(genes_table) %>%
  select(locus_tag,id)

locus_tags_aliases = read_csv('data/harbison/locus_tags_aliases.csv') %>%
  left_join(genes_table, by = c('curr_notation'='locus_tag')) %>%
  dplyr::rename(locus_tag=harbison_locus_tag) %>%
  select(locus_tag,id)

locus_tags_complete = locus_tags %>%
  filter(!is.na(id)) %>%
  bind_rows(locus_tags_aliases)

combined_harbison_with_id = combined_harbison %>%
  left_join(locus_tags_complete) %>%
  rename(gene_id=id)

tf_name_df = combined_harbison %>%
  select(tf) %>%
  mutate(tf = trimws(toupper(tf))) %>%
  distinct() %>%
  left_join(genes_table, by = c('tf' = 'gene')) %>%
  select(tf,id)

locus_tag_tfs = tf_name_df %>%
  filter(!complete.cases(.)) %>%
  select(tf) %>%
  left_join(genes_table, by = c('tf' = 'locus_tag')) %>%
  filter(!is.na(id)) %>%
  select(tf,id)

tf_name_df = tf_name_df %>%
  filter(!tf %in% locus_tag_tfs$tf) %>%
  bind_rows(locus_tag_tfs)

tf_name_df_na = read_csv("data/harbison/tf_name_aliases.csv") %>%
  left_join(genes_table) %>%
  select(tf,id)

tf_name_df = tf_name_df %>%
  filter(complete.cases(.)) %>%
  bind_rows(tf_name_df_na)

combined_harbison_with_id_tf_id = combined_harbison_with_id %>%
  filter(tf != 'A1 (MATA1)') %>%
  mutate(tf = trimws(toupper(tf))) %>%
  left_join(tf_name_df) %>%
  dplyr::rename(tf_id=id)

# combined_harbison_with_id_tf_id %>%
#   select(tf_id,cond,gene_id,binding_ratio,pval) %>%
#   group_by(tf_id,cond) %>%
#   group_walk(~{
#     name=paste(.y[1,], collapse='_')
#     message(name)
#     write_csv(.x,file.path(here('data/harbison/by_tf'),
#                            paste0(name,'.csv.gz')))
#   })



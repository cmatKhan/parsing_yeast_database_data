library(tidyverse)
library(here)

genes_table = read_csv(here('data/genes_table.csv.gz'))

hu_data = list(
  de = here("data/hu/gkq232_SuppTable1C_KOTargetGenes_Matrix_LogFoldChange.dat.xz"),
  pval = here("data/hu/gkq232_SuppTable1B_KOTargetGenes_Matrix_PValue.dat.xz")
)

hu_de_target_genes = as.character(read_delim(hu_data$de,
                             delim=' ',
                             col_names = FALSE,
                             n_max = 1)[1,])

hu_de_colnames = c('tf', hu_de_target_genes)

de_df = read.delim(
  hu_data$de,
  sep=' ',
  col.names = hu_de_colnames,
  skip=1) %>%
  as_tibble()

hu_pval_target_genes = as.character(read_delim(hu_data$pval,
                                             delim=' ',
                                             col_names = FALSE,
                                             n_max = 1)[1,])

hu_pval_colnames = c('tf', hu_pval_target_genes)

pval_df = read.delim(
  hu_data$pval,
  sep=' ',
  col.names = hu_pval_colnames,
  skip=1) %>%
  as_tibble()

hu_df = de_df %>%
  pivot_longer(-tf,names_to='target_gene',values_to='effect') %>%
  left_join(pval_df %>%
              pivot_longer(-tf,names_to='target_gene',values_to='pval')) %>%
  mutate(tf = toupper(tf),
         target_gene = str_replace(target_gene,'\\.','-')) %>%
  # remove deleted orfs
  filter(!target_gene %in% c('YER187W-A', 'YCR103C', 'YCL006C'))

gene_name_tf_map = hu_df %>%
  distinct(tf) %>%
  left_join(genes_table %>%
              select(gene,id) %>%
              dplyr::rename(tf_id=id),
            by = c('tf' = 'gene'))

locus_tag_tf_map = gene_name_tf_map %>%
  filter(is.na(tf_id)) %>%
  select(-tf_id) %>%
  left_join(genes_table %>%
              select(locus_tag,id) %>%
              dplyr::rename(tf_id=id),
            by = c('tf' = 'locus_tag'))

alias_tf_map = tibble(
  tf =    c('CAF17','RCS1','RLR1','ZMS1','RIS1'),
  tf_id = c(3839,    2419, 5658,  3844,   6310)
)

tf_map = bind_rows(
  gene_name_tf_map,
  locus_tag_tf_map,
  alias_tf_map
) %>%
  filter(!is.na(tf_id))

gene_map = hu_df %>%
  distinct(target_gene) %>%
  left_join(select(genes_table, locus_tag, id)
            %>% dplyr::rename(gene_id=id),
            by = c('target_gene' = 'locus_tag'))

alias_gene_map = gene_map %>%
  filter(is.na(gene_id)) %>%
  select(-gene_id) %>%
  left_join(select(genes_table, alias, id)
            %>% dplyr::rename(gene_id=id),
            by = c('target_gene' = 'alias')) %>%
  filter(!is.na(gene_id))

by_hand_gene_map = read_csv('data/hu/deleted_merged_gene_map.csv') %>%
  left_join(select(genes_table, locus_tag,id) %>% dplyr::rename(gene_id=id)) %>%
  select(target_gene,gene_id)

final_gene_map = gene_map %>%
  filter(complete.cases(.)) %>%
  bind_rows(alias_gene_map %>% filter(complete.cases(.))) %>%
  bind_rows(by_hand_gene_map)

setdiff(unique(hu_df$target_gene), final_gene_map$target_gene)

hu_df_with_ids = hu_df %>%
  left_join(tf_map) %>%
  left_join(final_gene_map) %>%
  select(tf_id,gene_id,effect,pval) %>%
  group_by(tf_id)

# hu_df_with_ids %>%
#   group_walk(~write_csv(.x,file.path(here('data/hu/tf_split'),
#                                     paste0(.y,'.csv.gz'))))
#

library(tidyverse)
library(here)

gene_table = read_csv("data/genes_table.csv")

# see the README in the mcisaac data directory
df = read_tsv(here("data/mcisaac/idea_tall_expression_data.tsv")) %>%
  filter(!TF %in% c("GEV","Z3EV")) %>%
  mutate(GeneName = str_remove(GeneName,',')) %>%
  mutate(GeneName = str_remove(GeneName, '\''))

tf_table = df %>%
  select(TF) %>%
  distinct() %>%
  filter(TF != 'YLL054C') %>%
  left_join(gene_table, by = c('TF' = 'gene')) %>%
  select(TF,id) %>%
  bind_rows(tibble(TF='YLL054C', id=4270)) %>%
  dplyr::rename(tf_id=id)

tmp = df %>%
  select(GeneName) %>%
  distinct() %>%
  left_join(gene_table, by = c('GeneName' = 'gene')) %>%
  filter(complete.cases(.)) %>%
  bind_rows(df %>%
              select(GeneName) %>%
              distinct() %>%
              left_join(gene_table, by = c('GeneName' = 'locus_tag')) %>%
              filter(complete.cases(.))) %>%
  select(GeneName, id) %>%
  dplyr::rename(gene_id=id)

# df %>%
#   select(GeneName) %>%
#   distinct(GeneName) %>%
#   filter(!GeneName %in% tmp$GeneName) %>%
#   mutate(id=NA) %>%
#   select(GeneName, id) %>%
#   dplyr::rename(gene_id=id) %>%
#   write_csv("data/mcisaac/mcisaac_gene_lookup.txt")

mcisaac_gene_aliases = read_csv("data/mcisaac/mcisaac_gene_lookup.txt") %>%
  left_join(gene_table) %>%
  dplyr::rename(gene_id=id) %>%
  select(GeneName,gene_id)

mcisaac_gene_table = tmp %>%
  bind_rows(mcisaac_gene_aliases)

mcisaac_with_ids = df %>%
  left_join(tf_table) %>%
  left_join(mcisaac_gene_table)

x = mcisaac_with_ids %>%
  select(tf_id,gene_id,strain,date,restriction,mechanism,time,starts_with('log2')) %>%
  group_by(tf_id,strain,date,restriction,mechanism,time) %>%
  group_split() %>%
  .[[1]]

library(tidyverse)
library(httr)
library(here)

source(here('R/get_combined_data.R'))

url_generator = function(endpiont){
  base = "http://127.0.0.1:8000/api/"
  paste0(base, endpiont)
}

token = Sys.getenv('local_token')

header <- httr::add_headers(
  "Authorization" = paste("token", token, sep = " ")
)

genomicfeature_df = GET(file.path(url_generator('genomicfeature'),'export/'), header) %>%
  httr::content() %>%
  select(-c(uploader_id,upload_date,modifier_id,modified_date))

regulator_df = GET(file.path(url_generator('regulator'),'export/'), header) %>%
  httr::content() %>%
  select(-c(uploader_id,upload_date,modifier_id,modified_date)) %>%
  left_join(genomicfeature_df, by = c('genomicfeature_id' = 'id'))

fileformat_df = GET(file.path(url_generator('fileformat'),'export/'), header) %>%
  httr::content() %>%
  select(-c(uploader_id,upload_date,modifier_id,modified_date))

datasource_df = GET(file.path(url_generator('datasource'),'export/'), header) %>%
  httr::content() %>%
  select(-c(uploader_id,upload_date,modifier_id,modified_date)) %>%
  left_join(fileformat_df, by = c('fileformat_id' = 'id'))


binding_df = GET(file.path(url_generator('binding'),'export/'), header) %>%
  httr::content() %>%
  select(-c(uploader_id,upload_date,modifier_id,modified_date)) %>%
  left_join(datasource_df, by = c('source_id' = 'id')) %>%
  left_join(regulator_df, by = c('regulator_id' = 'id'))


expression_df = GET(file.path(url_generator('expression'),'export/'), header) %>%
  httr::content() %>%
  select(-c(uploader_id,upload_date,modifier_id,modified_date)) %>%
  left_join(datasource_df, by = c('source_id' = 'id')) %>%
  left_join(regulator_df, by = c('regulator_id' = 'id')) %>%
  mutate(file = file.path('/home/oguzkhan/code/yeastregulatorydb/yeastregulatorydb/media', file))

promotersetsig_df = GET(file.path(url_generator('promotersetsig'),'export/'), header) %>%
  httr::content() %>%
  select(-c(uploader_id,upload_date,modifier_id,modified_date)) %>%
  left_join(fileformat_df, by = c('fileformat_id' = 'id')) %>%
  left_join(binding_df %>% select(id, name, regulator_symbol, condition, replicate, batch), by = c('binding_id' = 'id')) %>%
  mutate(file = file.path('/home/oguzkhan/code/yeastregulatorydb/yeastregulatorydb/media', file))

# number of regulators per datasource: binding
# ------------------------------
binding_df %>%
  group_by(regulator_symbol, name) %>%
  tally() %>%
  group_by(name) %>%
  tally()


# number of regulators per datasource: expression
# ------------------------------
expression_df %>%
  group_by(regulator_symbol, name) %>%
  tally() %>%
  group_by(name) %>%
  tally()

# Overlap btwn binding and expression by datasource
# ------------------------------
expression_df %>%
  group_by(regulator_symbol, name) %>%
  tally() %>%
  select(-n) %>%
  dplyr::rename(expression_source = name) %>%
  left_join(
    binding_df %>%
      group_by(regulator_symbol, name) %>%
      tally() %>%
      select(-n) %>%
      dplyr::rename(binding_source = name),
    relationship = 'many-to-many'
  ) %>%
  ungroup() %>%
  filter(complete.cases(.)) %>%
  group_by(expression_source,binding_source) %>%
  tally()

expression_tfs = expression_df %>%
  group_by(regulator_symbol, name) %>%
  tally() %>%
  select(-n) %>%
  dplyr::rename(expression_source = name) %>%
  left_join(
    binding_df %>%
      group_by(regulator_symbol, name) %>%
      tally() %>%
      select(-n) %>%
      dplyr::rename(binding_source = name),
    relationship = 'many-to-many'
  ) %>%
  ungroup() %>%
  filter(complete.cases(.)) %>%
  group_by(expression_source,binding_source) %>%
  ungroup() %>%
  distinct(regulator_symbol) %>%
  pull()

# top30_list = map(expression_tfs, summarize_binding_response_tf)
#
# top30_df = bind_rows(top30_list)
#
# top30_df = top30_df %>%
#   select(-ends_with('.y'))
# colnames(top30_df) = str_remove(colnames(top30_df), '.x$')
# write_csv(top30_df,"data/test_data/top_30_all_expr_tfs_20230119.csv.gz")
top30_df = read_csv("data/test_data/top_30_all_expr_tfs_20230119.csv.gz")

# x %>%
#   group_by(regulator_symbol, promotersetsig_id, binding_source, expression_source, responsive) %>%
#   tally() %>%
#   mutate(y = sum(n), ratio = n / sum(n),
#          binding_assay = case_when(
#            startsWith(binding_source, 'brent') ~ 'nf_callingcards',
#            startsWith(binding_source, "mitra") ~ "mitra_callingcards",
#            startsWith(binding_source, 'chipexo') ~ 'chipexo',
#            startsWith(binding_source, 'harb') ~ 'harbison',
#            .default = 'unknown'
#          )) %>%
#   filter(responsive) %>%
#   ggplot(aes(expression_source, ratio, fill=binding_assay, color = promotersetsig_id)) +
#   geom_col(position='dodge') +
#   facet_grid(~regulator_symbol) +
#   theme_bw() +
#   theme(text = element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1)) +
#   ggtitle('Fraction of the top 30 genes by binding pvalue that are responsive') +
#   guides(color = 'none')

top30_df %>%
  mutate(promotersetsig_id = factor(promotersetsig_id)) %>%
  ungroup() %>%
  filter(complete.cases(.)) %>%
  group_by(regulator_symbol, promotersetsig_id, binding_source, expression_source) %>%
  arrange(binding_pvalue) %>%
  slice_head(n=30) %>%
  group_by(regulator_symbol, promotersetsig_id, binding_source, expression_source, responsive) %>%
  tally() %>%
  mutate(y = sum(n), ratio = n / sum(n),
         binding_assay = case_when(
           startsWith(binding_source, 'brent') ~ 'nf_callingcards',
           startsWith(binding_source, "mitra") ~ "mitra_callingcards",
           startsWith(binding_source, 'chipexo') ~ 'chipexo',
           startsWith(binding_source, 'harb') ~ 'harbison',
           .default = 'unknown'
         )) %>%
  filter(responsive) %>%
  ggplot(aes(expression_source, ratio, fill = binding_assay)) +
  geom_boxplot() +
  theme_bw() +
  theme(text = element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('Distributions of responsive ratio across TFs',
          subtitle = 'top 30 genes by binding pval') +
  guides(color = 'none')


go_res_df = x %>%
  ungroup() %>%
  filter(complete.cases(.)) %>%
  group_by(regulator_symbol, promotersetsig_id, expression_source) %>%
  arrange(regulator_symbol, promotersetsig_id,expression_source,binding_pval) %>%
  left_join(select(genes_table,id,locus_tag), by = c('gene_id' = 'id'))

go_results = go_res_df %>%
  group_map( ~ {message('here!')
                gprofiler2::gost(.$locus_tag,
                    organism='scerevisiae',
                    ordered_query=TRUE)})

names(go_results) = paste(group_keys(go_res_df)$binding_source,
                          group_keys(go_res_df)$expression_source, sep = ';')

go_res_df_list = map(go_results, ~{
  df = .$result %>%
    group_by(source) %>%
    arrange(p_value) %>%
    slice_head(n=1)}) %>%
  bind_rows(.id = 'tmp') %>%
  separate(tmp, c('binding_source','expression_source'), sep = ';')

go_res_df_list %>%
  filter(term_size >= 3 & term_size <= 200 ) %>%
  mutate(binding_assay = case_when(
    startsWith(binding_source, 'cc') ~ 'callingcards',
    startsWith(binding_source, 'chipexo') ~ 'chipexo',
    startsWith(binding_source, 'harb') ~ 'harbison',
    .default = 'unknown'
  )) %>%
  ungroup() %>%
  group_by(binding_source,expression_source,source) %>%
  slice_max(p_value, n=1) %>%
  ggplot(aes(expression_source, -log10(p_value), fill = binding_assay, color = binding_source)) +
  geom_col(position='dodge') +
  theme_bw() +
  facet_grid(~source) +
  theme(text = element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('GSEA by binding pvalue')


go_results_ora = go_res_df %>%
  ungroup() %>%
  group_by(binding_source, expression_source) %>%
  slice_max(-binding_pval, n=500) %>%
  group_map( ~ gost(.$locus_tag,
                    organism='scerevisiae'))

names(go_results_ora) = paste(group_keys(go_res_df %>% ungroup() %>% group_by(binding_source,expression_source))$binding_source,
                              group_keys(go_res_df %>% ungroup() %>% group_by(binding_source,expression_source))$expression_source, sep = ';')

go_res_df_list_ora = map(go_results_ora, ~{
  df = .$result %>%
    group_by(source) %>%
    arrange(p_value) %>%
    slice_head(n=1)}) %>%
  bind_rows(.id = 'tmp') %>%
  separate(tmp, c('binding_source','expression_source'), sep = ';')

go_res_df_list_ora %>%
  filter(term_size >= 3 & term_size <= 200 ) %>%
  mutate(binding_assay = case_when(
    startsWith(binding_source, 'cc') ~ 'callingcards',
    startsWith(binding_source, 'chipexo') ~ 'chipexo',
    startsWith(binding_source, 'harb') ~ 'harbison',
    .default = 'unknown'
  )) %>%
  ungroup() %>%
  group_by(binding_source,expression_source,source) %>%
  slice_max(p_value, n=1) %>%
  ggplot(aes(expression_source, -log10(p_value), fill = binding_assay, color = binding_source)) +
  geom_col(position='dodge') +
  theme_bw() +
  facet_grid(~source) +
  theme(text = element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('ORA of top 50 genes by binding pvalue')


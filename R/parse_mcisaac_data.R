library(tidyverse)
library(here)
library(httr)
library(jsonlite)

gene_table = read_csv("data/genome_files/genes_table.csv.gz")

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
  mutate(date=str_remove_all(as.Date(date, format='%m/%d/%Y'),'-')) %>%
  left_join(tf_table) %>%
  left_join(mcisaac_gene_table) %>%
  select(TF, tf_id,gene_id,strain,date,restriction,
         mechanism,time,starts_with('log2'))


x = mcisaac_with_ids %>%
  group_by(tf_id,strain,date,restriction,mechanism,time) %>%
  group_split() %>%
  .[[1]]

# mcisaac_with_ids %>%
#   mutate(
#     across(starts_with("log2_"), ~format(., nsmall = 1))
#   ) %>%
#   group_by(tf_id,strain,date,restriction,mechanism,time) %>%
#   group_walk(~{
#     write_csv(.x,
#               file.path('data/mcisaac/by_tf',
#                            paste0(paste(.y[1,],collapse='_'),
#                                   '.csv.gz')))})

mcisaac_with_ids_upload_df = mcisaac_with_ids %>%
  group_by(tf_id, TF, strain, date, restriction, mechanism, time) %>%
  tally() %>%
  mutate(file = file.path('data/mcisaac/by_tf',
                          paste0(tf_id, "_", strain, "_", date, "_",
                                 restriction, "_", mechanism, "_", time, ".csv.gz")))


prepare_http_upload <- function(row) {

  stopifnot(file.exists(row$file))

  # Extract data from the row
  list(
    mechanism = tolower(row$mechanism),
    restriction = row$restriction,
    time = row$time,
    file = httr::upload_file(row$file,type="application/gzip"),
    regulator_symbol = row$TF,
    strain = row$strain,
    source_name = 'mcisaac_oe',
    notes = paste0('strain_id:', row$strain, ' ; ', 'date:', row$date)
  )


}

url = "http://127.0.0.1:8000/api/expression/"

token = Sys.getenv('local_token')

header <- httr::add_headers(
"Authorization" = paste("token", token, sep = " ")
)

# response_list = map(1:nrow(mcisaac_with_ids_upload_df), function(i) {
#   message(i)
#   body_data = prepare_http_upload(mcisaac_with_ids_upload_df[i,])
#   res = POST(url, config = timeout(180), header, body = body_data, encode = "multipart")
#   message(httr::status_code(res))
#   res
# })


library(tidyverse)
library(here)

source(here('R/post_data.R'))

# note that the token is saved in token_local
# response = POST('http://127.0.0.1:8000/auth-token/',
#                 body = list(username='username',password='password'))

token = Sys.getenv('local_token')

chrmap_raw = read_csv(here("data/chrmap.csv.gz"))

chrmap_upload = chrmap_raw %>%
  select(-c(uploader,uploadDate,modifiedBy,modified,id)) %>%
  filter(!numbered %in% c(1,2, 3, 4, 5, 6, 7, 8))

chr_url = 'http://127.0.0.1:8000/api/chrmap/'

#post_data(chrmap_upload, chr_url, token)

gene_raw = read_csv("data/genes_table.csv.gz")

gene_upload = gene_raw %>%
  dplyr::rename(biotype=gene_biotype, symbol=gene) %>%
  select(-c(uploader,uploadDate,modifiedBy,modified, biotype)) %>%
  mutate(start=ifelse(id==7150,1,start),
         end=ifelse(id==7150,2,end)) %>%
  filter(!locus_tag %in% c('YAL069W'))

genomic_feature_url = "http://127.0.0.1:8000/api/genomicfeature/"

post_data(gene_upload, genomic_feature_url, token)

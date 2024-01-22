library(tidyverse)
library(here)

# NOTE!!  the associated_feature must be to the genomicfeature `id`

promoter_regions = read_csv(here("data/promoter_regions.csv.gz"))

yiming_ucsc = promoter_regions %>%
  filter(source=='yiming') %>%
  select(chr,start,end,associated_feature,score,strand) %>%
  dplyr::rename(name=associated_feature) %>%
  left_join(chrmap %>% select(id,ucsc,seqlength), by =c('chr'='id')) %>%
  mutate(chr=ucsc) %>%
  mutate(end = ifelse(end>seqlength, seqlength,end)) %>%
  select(-c(ucsc,seqlength))

write_tsv(yiming_ucsc, here("data/promoters/yiming_ucsc.bed.gz"))

# note -- checked by hand. this does not have
# any endpoints that exceed the seqlengths
not_orf_ucsc = promoter_regions %>%
  filter(source=='not_orf') %>%
  select(chr,start,end,associated_feature,score,strand) %>%
  dplyr::rename(name=associated_feature) %>%
  left_join(chrmap %>% select(id,ucsc,seqlength), by =c('chr'='id')) %>%
  mutate(chr=ucsc) %>%
  select(-ucsc,seqlength)

write_tsv(not_orf_ucsc, here("data/promoters/not_orf_ucsc.bed.gz"))

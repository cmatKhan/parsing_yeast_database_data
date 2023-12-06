library(tidyverse)
library(here)

chr_map = read_csv('data/chrmap.csv.gz')

chipexo = read_csv(here("data/chip_exo/parsed/10352.csv.gz"))

promoter_regions = read_csv("data/promoter_regions.csv.gz")  %>%
  filter(source=='yiming') %>%
  left_join(chr_map, by = c('chr'='id')) %>%
  select(numbered,start,end,associated_feature,score,strand) %>%
  dplyr::rename(chr=numbered) %>%
  mutate(chr=paste0('chr',chr))

promoter_regions %>%
  inner_join(chipexo, by = 'chr')  %>%
  filter(coord >= start & coord <= end) %>%
  group_by(chr, start, end, associated_feature)


chipexo_test = tibble(
  chr=rep('chr1',3),
  coord=c(2190, 2200, 2500),
  YPD_Log2Fold=c(3,6, 0),
  YPD_Log2P=c(-300,-600, 0)
)

promoter_regions_test = promoter_regions %>%
  filter(associated_feature %in% c(3,4)) %>%
  inner_join(chipexo_test, relationship = 'many-to-many')%>%
  filter(coord >= start & coord <= end) %>%
  group_by(chr, start, end, associated_feature) %>%
  summarize(n_sig_peaks = n(),
            max_fc = max(YPD_Log2Fold),
            min = min(YPD_Log2P),
            .groups = 'drop_last')

sig_promoter_regions = promoter_regions %>%
  inner_join(chipexo, relationship = 'many-to-many') %>%
  filter(coord >= start & coord <= end) %>%
  group_by(chr, start, end, associated_feature) %>%
  summarize(n_sig_peaks = n(),
            max_fc = max(YPD_log2Fold),
            min_pval = min(YPD_log2P),
            .groups = 'drop_last')

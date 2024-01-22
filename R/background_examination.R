library(tidyverse)
library(here)

yiming_promoters = read_tsv(here("data/promoters/yiming_ucsc.bed.gz"))
adh1_background  = read_tsv("~/code/yeastregulatorydb/yeastregulatorydb/regulatory_data/tests/test_data/background/adh1_background.qbed.gz") %>%
  distinct(chr,start,end)
hap5_17 = read_csv(here("data/hap5_testset/cc_expr17.qbed_sig"))

promoter_background_joined = yiming_promoters %>%
  dplyr::rename(promoter_start = start,promoter_end = end) %>%
  mutate(promoter_width = promoter_end-promoter_start) %>%
  left_join(adh1_background, by = "chr", relationship = 'many-to-many') %>%
  filter(promoter_start <= start & promoter_end >= end)

background_promoter_analysis_df = promoter_background_joined %>%
  group_by(chr,promoter_start,promoter_end, strand) %>%
  reframe(hops = n()) %>%
  mutate(promoter_occupancy = hops/(promoter_end-promoter_start))

background_promoter_analysis_df %>%
  ggplot(aes(x=promoter_occupancy)) +
  geom_histogram(bins = 100) +
  facet_wrap(~strand)

background_promoter_analysis_df %>%
  ggplot(aes(hops)) +
  geom_histogram(bins = 100) +
  facet_wrap(~strand)

hap5_17_promoter_analysis_df = hap5_17 %>%
  mutate(promoter_width = end-start,
         promoter_occupancy = experiment_hops/promoter_width)

hap5_17_promoter_analysis_df %>%
  ggplot(aes(x=promoter_occupancy)) +
  geom_histogram(bins = 100) +
  facet_wrap(~strand)

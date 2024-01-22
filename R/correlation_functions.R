library(tidyverse)
library(here)
library(corrplot)
library(gprofiler2)

#https://chat.openai.com/share/28b07237-540d-4bec-a23f-fa71f4394eef

chrmap =read_csv("data/genome_files/chrmap.csv.gz")
promoter_df = read_tsv("data/promoters/yiming_ucsc.bed.gz")
genes_df = read_csv("data/genome_files/genes_table.csv.gz")

translate_chr_format = function(df){
  df %>%
    left_join(chrmap %>% select(id,ucsc, seqlength), by =c('chr'='id')) %>%
    mutate(chr=ucsc) %>%
    mutate(end = ifelse(end>seqlength, seqlength,end)) %>%
    select(-c(ucsc,seqlength))
}

files = list.files(here("data/hap5_testset/"), full.names = TRUE)

cc_data = map(files[str_detect(files,'.qbed_sig')], read_csv)

names(cc_data) = basename(files[str_detect(files,'.qbed_sig')]) %>% str_remove('.qbed_sig')

cc_data_df = cc_data %>%
  bind_rows(.id='ccexperiment') %>%
  dplyr::rename(gene_id = name)

chipexo_data = map(files[str_detect(files,'csv.gz_sig')], read_csv)

names(chipexo_data) = basename(files[str_detect(files,'csv.gz_sig')]) %>% str_remove('.csv.gz_sig')

chipexo_data_df = chipexo_data %>%
  bind_rows(.id='chipexoexperiment') %>%
  dplyr::rename(gene_id = name)

harbison = read_csv("data/hap5_testset/harbison_6498.csv.gz") %>%
  mutate(source = 'harbison') %>%
  group_by(gene_id, source) %>%
  summarize(effect = max(effect), pval = min(pval)) %>%
  ungroup()

binding_df = harbison %>%
  bind_rows(
    cc_data_df %>%
      select(gene_id,ccexperiment,callingcards_enrichment,poisson_pval) %>%
      dplyr::rename(effect = callingcards_enrichment,pval = poisson_pval,source = ccexperiment)
  ) %>%
  bind_rows(
    chipexo_data_df %>%
      select(gene_id,chipexoexperiment,max_fc,min_pval) %>%
      dplyr::rename(effect = max_fc, pval = min_pval,source = chipexoexperiment)
  )

expression_data = map(files[str_detect(files,'kemmeren|hu|mcisaac')], read_csv)

names(expression_data) = basename(files[str_detect(files,'kemmeren|hu|mcisaac')]) %>% str_remove('.csv.gz_sig')

expression_df = expression_data$mcisaac_6498_SMY2076_20160321_P_ZEV_15.csv.gz %>%
  select(gene_id, log2_shrunken_timecourses) %>%
  dplyr::rename(effect = log2_shrunken_timecourses) %>%
  mutate(pval=0, source='mcisaac') %>%
  group_by(gene_id,source) %>%
  summarize(effect = max(effect), pval = min(pval)) %>%
  ungroup() %>%
  bind_rows(
    expression_data$kemmeren_6498.csv.gz %>%
      select(gene_id, Madj, pval) %>%
      dplyr::rename(effect=Madj) %>%
      mutate(source = 'kemmeren')
  ) %>%
  bind_rows(expression_data$hu_6498.csv.gz %>%
              mutate(source = 'hu') %>%
              group_by(gene_id,source) %>%
              summarize(effect = max(effect), pval = min(pval)))

full_df = binding_df %>%
  bind_rows(expression_df)

effect_df = full_df %>%
  select(gene_id,effect,source) %>%
  pivot_wider(id_cols=gene_id, names_from = source, values_from = effect)




data = effect_df %>%
  select(-chipexo_hap5_rep1_10535, -chipexo_hap5_rep2_10706) %>%
  filter(complete.cases(.)) %>%
  select(-gene_id)

# Assuming you have a dataframe 'data' with each column representing an assay
pca_result <- prcomp(data, center = TRUE, scale. = TRUE)

# Plotting the first two principal components
plot(pca_result$x[,1], pca_result$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA of Assay Data")
text(pca_result$x[,1], pca_result$x[,2], labels = rownames(pca_result$x), pos = 4)

# Calculating variance explained by each principal component
var_explained <- pca_result$sdev^2 / sum(pca_result$sdev^2)

# Creating a scree plot
plot(var_explained, xlab = "Principal Component", ylab = "Proportion of Variance Explained",
     type = "b", pch = 19, main = "Scree Plot")

# Creating a biplot
biplot(pca_result, scale = 0, cex = 0.6, main = "PCA Biplot")


## corrplot CC data replicates

correlation_data = cc_data_df %>%
  select(ccexperiment, gene_id, poisson_pval) %>%
  group_by(ccexperiment) %>%
  arrange(poisson_pval) %>%
  mutate(rank_raw_order = row_number(),
         rank_func = log(rank(-poisson_pval, ties.method="average")))

correlation_data %>%
  select(ccexperiment, gene_id, poisson_pval) %>%
  pivot_wider(id_cols=gene_id, names_from = ccexperiment, values_from = poisson_pval) %>%
  select(-gene_id) %>%
  cor(use='pairwise.complete.obs') %>%
  corrplot(method='number', type='upper', order='AOE', title='correlation by pvalue')

correlation_data %>%
  pivot_wider(id_cols=gene_id, names_from = ccexperiment, values_from = rank_func) %>%
  select(-gene_id) %>%
  cor(use='pairwise.complete.obs') %>%
  corrplot(method='number', type='upper', order='AOE', title='correlation by rank')

# Assuming you have a dataframe 'data' with each column representing an assay

pca_result_pval = correlation_data %>%
  filter(poisson_pval < 0.05) %>%
  pivot_wider(id_cols=gene_id, names_from = ccexperiment, values_from = poisson_pval) %>%
  select(-gene_id) %>%
  filter(complete.cases(.)) %>%
  prcomp(center = TRUE, scale. = TRUE)

# Plotting the first two principal components
plot(pca_result_pval$x[,1], pca_result_pval$x[,2], xlab = "PC1", ylab = "PC2", main = "PCA of Assay Data")
text(pca_result_pval$x[,1], pca_result_pval$x[,2], labels = rownames(pca_result_pval$x), pos = 4)

# Calculating variance explained by each principal component
var_explained <- pca_result_pval$sdev^2 / sum(pca_result_pval$sdev^2)

# Creating a scree plot
plot(var_explained, xlab = "Principal Component", ylab = "Proportion of Variance Explained",
     type = "b", pch = 19, main = "Scree Plot")

# Creating a biplot
biplot(pca_result_pval, scale = 0, cex = 0.6, main = "PCA Biplot")

# rank analysis
#
# We want some things that we can represent as a single number, such as the
# fraction of the top N that are responsive, where N is a user-changeable
# parameter defaulting to maybe 30.  (3) We pick a way of displaying them.
# E.g. for all binding versus all response for a single TF, we could display
# them as a clustered bar chart (clustered by response dataset, colored by
# binding dataset). (4) We make those displays for a bunch of TFs and decide
# whether we like them. (edited)

x = binding_df %>%
  dplyr::rename(binding_source = source,
                binding_effect=effect,
                binding_pval=pval) %>%
  left_join(expression_df %>%
              dplyr::rename(expression_source = source,
                            expression_effect = effect,
                            expression_pval = pval), by =c('gene_id'),
            relationship = 'many-to-many') %>%
  mutate(responsive = case_when(
    expression_source == 'mcisaac' & expression_effect > 0 ~ TRUE,
    expression_source == 'kemmeren' & expression_pval < 0.5  ~ TRUE,
    expression_source == 'hu' & expression_pval < 0.5 ~ TRUE,
    .default = FALSE
  ))

x %>%
  ungroup() %>%
  filter(complete.cases(.)) %>%
  group_by(binding_source, expression_source) %>%
  top_n(-30, binding_pval) %>%
  group_by(binding_source, expression_source, responsive) %>%
  tally() %>%
  mutate(y = sum(n), ratio = n / sum(n),
         binding_assay = case_when(
           startsWith(binding_source, 'cc') ~ 'callingcards',
           startsWith(binding_source, 'chipexo') ~ 'chipexo',
           startsWith(binding_source, 'harb') ~ 'harbison',
           .default = 'unknown'
         )) %>%
  filter(responsive) %>%
  ggplot(aes(expression_source, ratio, fill=binding_assay, color = binding_source)) +
  geom_col(position='dodge') +
  theme_bw() +
  theme(text = element_text(size=18), axis.text.x = element_text(angle = 45, hjust = 1)) +
  ggtitle('Fraction of the top 30 genes by binding pvalue that are responsive')


# Next, we go back to (2) and repeat. I think it would be nice to have some
# generic ways of comparing scores that don’t depend on having one binding and
# one response dataset but would work on any pair of datasets. One obvious way
# of doing that would be using correlations. Here I think we might need to try
# a few different methods, such as Pearson correlation of raw scores and
# log-transformed scores (oriented so larger scores indicate more binding or
# response) and rank correlation using untransformed or log transformed ranks.
# Then we go back to (3), pick some way of displaying them. Here I think, in
# addition to TF-by-TF comparisons, we should also look at some global
# dataset-by-dataset displays. Easiest is probably box and whisker displaying
# the distribution for each pair of datasets. (4) Look at these and decide if
# we like them and whether we want to keep all as options with one as default,
# or just pick one.

## GO metrics

go_res_df = x %>%
  ungroup() %>%
  filter(complete.cases(.)) %>%
  group_by(binding_source, expression_source) %>%
  arrange(binding_source,expression_source,binding_pval) %>%
  left_join(select(genes_df,id,locus_tag), by = c('gene_id' = 'id'))

go_results = go_res_df %>%
  group_map( ~ gost(.$locus_tag,
                    organism='scerevisiae',
                    ordered_query=TRUE))

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

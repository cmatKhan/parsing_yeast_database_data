library(tidyverse)
library(here)

add_datatype_to_colnames = function(df, skip_indicies){
  # Suffixes to append
  suffixes <- c("_M", "_A", "_pval")

  # Repeat the suffixes to match the length of my_vector
  repeated_suffixes <- rep(suffixes, length.out = length(colnames(df)[-skip_indicies]))

  # Append the suffixes to each element of my_vector
  modified_vector <- paste0(colnames(df)[-skip_indicies], repeated_suffixes)

  colnames(df)[-skip_indicies] = modified_vector

  df[-1,]
}

# note that these are created after the `delete_all_mutants_controls`
# dataframe. I created that df, used it to create the tf and gene maps
# saved them and then use them to merge below
tf_map = read_csv(here('data/kemmeren/tf_map.csv.gz'))
# NOTE: this has a column added, wt_variable, which is derived from the
#wt_variable data from the kemmeren source data. This column can be used
#to remove wt_variable genes from the data with filter(!wt_variable)
gene_map = read_csv(here('data/kemmeren/gene_map.csv.gz'))
replicate_map = read_csv(here('data/kemmeren/replicate_table.csv.gz'))


# the headers vector was created with the code below from the raw data
# directly from the kemmeren site. It was used to read the full data set in
# and turn it into a parquet directory b/c the data itself was too large
# to store in github
deleteome_all_mutants_controls_headers =
  read_tsv(here('data/kemmeren/deleteome_all_mutants_controls.txt.xz'),
           n_max=1, name_repair='minimal') %>%
  add_datatype_to_colnames(skip_indicies=1:3) %>%
  colnames()

deleteome_all_mutants_controls_headers = str_replace(deleteome_all_mutants_controls_headers,
                                                     ' vs\\.? ', ';')

# this is necessary b/c there is one TF -- LUG1 and YLR352W -- which
# are actually the same thing (LUG1 is YLR352W). There are no other replicates.
# Regressing the Ms onto each other results in a significant association,
# but the effect is small (.08, pval .00132)
replicate_table = read_csv(here('data/kemmeren/replicate_table.csv.gz'))
# this takes minutes to execute and ~19GB
deleteome_all_mutants_controls =
  read.delim(here('data/kemmeren/deleteome_all_mutants_controls.txt.xz'),
           sep='\t',
           skip=2,
           check.names=FALSE,
           col.names=deleteome_all_mutants_controls_headers) %>%
  as_tibble() %>%
  pivot_longer(-c(reporterId, systematicName, geneSymbol),
               names_to='sample_metric', values_to='values') %>%
  separate(sample_metric, c('sample', 'metric'), sep="_") %>%
  pivot_wider(names_from='metric', values_from='values') %>%
  separate_wider_delim(cols=sample,
                       names=c('mutant', 'control'),
                       delim=";") %>%
  separate_wider_delim(cols=mutant,
                       names=c('tf', 'perturbation'),
                       delim="-",
                       too_few='error',
                       too_many='merge') %>%
  mutate(tf = toupper(str_replace(tf,',',''))) %>%
  filter(tf != 'WT') %>%
  left_join(tf_map) %>%
  left_join(gene_map) %>%
  left_join(replicate_table)

# replicate_table = deleteome_all_mutants_controls %>%
#   distinct(tf_id,tf) %>%
#   group_by(tf_id) %>%
#   mutate(replicate=row_number())
#
# by_hand_gene_map = tibble(
#   systematicName=c('YAR062W', 'YDL038C', 'snR10',     'YGR272C',   'YIL080W',    'YIL168W', 'YIR044C'),
#   locus_tag      =c('YAR061W', 'YDL039C', 'YNCG0013W', 'YGR271C-A', 'YIL082W-A', 'YIL167W', 'YIR043C')
# ) %>%
#   left_join(select(genes_table,id,locus_tag)) %>%
#   select(-locus_tag)
#
# tibble(
#   systematicName=unique(deleteome_all_mutants_controls$systematicName)) %>%
#   left_join(select(genes_table,id,locus_tag),
#             by=c('systematicName' = 'locus_tag')) %>%
#   filter(complete.cases(.)) %>%
#   bind_rows(by_hand_gene_map) %>%
#   dplyr::rename(gene_id=id) %>%
#   write_csv(here('data/kemmeren/gene_map.csv.gz'))

deleteome_all_mutants_svd_transformed_headers =
  read_tsv(here('data/kemmeren/deleteome_all_mutants_svd_transformed.txt.xz'),
           n_max=1,
           name_repair='minimal') %>%
  colnames()

deleteome_all_mutants_svd_transformed_headers[3:length(deleteome_all_mutants_svd_transformed_headers)] =
  paste0(deleteome_all_mutants_svd_transformed_headers[3:length(deleteome_all_mutants_svd_transformed_headers)],
         "_Madj")

deleteome_all_mutants_svd_transformed_headers =
  str_replace(deleteome_all_mutants_svd_transformed_headers,
              '\\.vs\\.\\.', ';')
deleteome_all_mutants_svd_transformed_headers =
  str_replace_all(deleteome_all_mutants_svd_transformed_headers,
              '\\.', '-')
deleteome_all_mutants_svd_transformed_headers[1] = 'systematicName'

deleteome_all_mutants_svd_transformed =
  read.delim(here('data/kemmeren/deleteome_all_mutants_svd_transformed.txt.xz'),
           sep='\t',
           check.names=FALSE,
           col.names=deleteome_all_mutants_svd_transformed_headers) %>%
  as_tibble() %>%
  # this locus is not in the `wt_control` (untransformed) data
  filter(commonName != 'Q0010') %>%
  # these two loci had the commonName entered as the systematicName
  mutate(systematicName =
           case_when(systematicName=='ANR2' ~ 'YKL047W',
                     systematicName=='CMS1' ~ 'YLR003C',
                     .default = systematicName)) %>%
  pivot_longer(-c(systematicName, commonName),
               names_to='sample_metric', values_to='values') %>%
  separate(sample_metric, c('sample', 'metric'), sep="_") %>%
  pivot_wider(names_from='metric', values_from='values') %>%
  separate_wider_delim(cols=sample,
                       names=c('mutant', 'control'),
                       delim=";") %>%
  separate_wider_delim(cols=mutant,
                       names=c('tf', 'perturbation'),
                       delim="-",
                       too_few='debug',
                       too_many='merge') %>%
  mutate(perturbuation=ifelse(!mutant_ok, FALSE, perturbation)) %>%
  select(-c(mutant_ok,mutant_pieces,mutant_remainder,mutant)) %>%
  mutate(tf = toupper(str_replace(tf,',',''))) %>%
  # the parsing incorrectly affects the `MF(ALPHA)1`, `MF(AlPHA)2` and
  # ARG5,6 loci
  mutate(tf = case_when(
    tf == 'MF' & perturbation == 'alpha-2-del' ~ 'MF(ALPHA)2',
    tf == 'MF' & perturbation == 'alpha-1-del' ~ 'MF(ALPHA)1',
    tf == 'ARG5' ~ 'ARG56',
    .default=tf)) %>%
  mutate(perturbation = case_when(
    perturbation == 'alpha-2-del' ~ 'del',
    perturbation == 'alpha-1-del' ~ 'del',
    perturbation == '6-del' ~ 'del',
    tf == 'YMR031C' ~ 'del',
    .default=perturbation)) %>%
  filter(tf != 'WT') %>%
  left_join(tf_map) %>%
  left_join(gene_map) %>%
  left_join(replicate_map)

x = deleteome_all_mutants_controls %>%
  left_join(select(
    deleteome_all_mutants_svd_transformed,systematicName,
    tf_id,gene_id,
    replicate,Madj)) %>%
  select(reporterId,tf_id,gene_id,replicate,perturbation,control,M,Madj,A,pval) %>%
  mutate(perturbation = 'del') %>%
  group_by(tf_id,gene_id,replicate,perturbation,control)%>%
  summarize(
    reporterId = reporterId[which.max(abs(M))],
    M = M[which.max(abs(M))],
    Madj = Madj[which.max(abs(Madj))],
    A = A[which.max(abs(A))],
    pval = min(pval),
    .groups = 'drop'
  ) %>%
  group_by(tf_id,replicate) %>%
  dplyr::rename(mechanism = perturbation) %>%
  mutate(mechanism = ifelse(mechanism=='del', 'TFKO', mechanism),
         control = str_replace(control,'-','_'))

# x %>%
#   group_walk(~{
#       message(paste(.y$tf_id,.y$replicate))
#       write_csv(select(.x,-c(perturbation,control)),file.path(here('data/kemmeren/by_tf'),
#                           paste0(paste(.y$tf_id,.y$replicate, sep="_"),
#                                  '.csv.gz')))})

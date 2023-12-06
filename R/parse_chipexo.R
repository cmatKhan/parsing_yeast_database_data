library(tidyverse)
library(here)

meta_df = read_tsv(here('data/chip_exo/SupplementaryData-Table4_Sample-Key_tabular.tab')) %>%
  dplyr::rename(chipexo_id = `Sample ID`,
                assay=`#Assay`,
                tf = `Yeast Target Common Name`,
                accession = Accession,
                sra_accession = `SRA Accession`,
                replicate=Replicate) %>%
  mutate(chipexo_id = as.character(chipexo_id),
         replicate = str_remove(replicate,'rep'))

genes_table = read_csv("data/genes_table.csv.gz")

chipexo_by_hand_gene_map = read_csv(here('data/chip_exo/chipexo_by_hand_gene_map.csv')) %>%
  left_join(genes_table) %>%
  select(tf, id)

chipexo_files = list.files(here('data/chip_exo/raw'), ".tabular$",
                           recursive = TRUE,
                           full.names = TRUE)

names(chipexo_files) = str_remove(basename(chipexo_files),
                                  '_chexmix_allevents.tabular')

read_in_chipexo = function(path){

  meta_line1=read_tsv(path, n_max = 1,skip=1)
  meta_line2=read_tsv(path, n_max = 1,skip=3)

  read_tsv(path,
           comment = '#',
           col_names = c('chr_coord', 'YPD_Sig','YPD_Ctrl','YPD_log2Fold','YPD_log2P','ActiveConds')) %>%
    separate(chr_coord, into = c('chr','coord'), remove=TRUE, sep=":") %>%
    mutate(sig_count = meta_line1$TotalSigCount,
           control_count = meta_line2$CtrlCount,
           sig_fraction = meta_line1$SignalFraction,
           sigCtrlScaling = meta_line2$SigCtrlScaling,
           cond = meta_line1$Name,
           parent_cond = meta_line2$ParentCond) %>%
    mutate(YPD_log2P = as.numeric(ifelse(YPD_log2P=='-Infinity', Inf, YPD_log2P)))

}

chipexo_df_list = map(chipexo_files, read_in_chipexo)

chipexo_df = chipexo_df_list %>%
  bind_rows(.id='chipexo_id') %>%
  left_join(meta_df) %>%
  mutate(assay = case_when(
    assay == "ChIP-exo" ~ 'chipexo',
    assay == "MNase-ChIP-seq" ~ 'mnase_chip_seq'),
    tf = toupper(tf)) %>%
  # CTDSER2 etc is modified RNA Pol II
  filter(assay == 'chipexo',
         !tf %in% c('CTDSER2', 'CTDSER5', 'CTDSER7')) %>%
  left_join(genes_table %>%
              select(id,gene) %>%
              dplyr::rename(tf_id=id),
            by = c('tf' = 'gene')) %>%
  left_join(chipexo_by_hand_gene_map) %>%
  mutate(tf_id = ifelse(is.na(tf_id), id, tf_id)) %>%
  select(-id)

# map2(chipexo_df_list, names(chipexo_df_list),
#      ~write_csv(select(.x,chr, coord, YPD_Sig, YPD_Ctrl, YPD_log2Fold, YPD_log2P),
#                 file.path('data/chip_exo/parsed', paste0(.y,'.csv.gz'))))



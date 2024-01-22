library(tidyverse)
library(here)
library(httr)

# meta_df = read_tsv(here('data/chip_exo/SupplementaryData-Table4_Sample-Key_tabular.tab')) %>%
#   dplyr::rename(chipexo_id = `Sample ID`,
#                 assay=`#Assay`,
#                 tf = `Yeast Target Common Name`,
#                 accession = Accession,
#                 sra_accession = `SRA Accession`,
#                 replicate=Replicate) %>%
#   mutate(chipexo_id = as.character(chipexo_id),
#          replicate = str_remove(replicate,'rep'))
#
# genes_table = read_csv("data/genome_files/genes_table.csv.gz")
#
# chipexo_by_hand_gene_map = read_csv(here('data/chip_exo/chipexo_by_hand_gene_map.csv')) %>%
#   left_join(genes_table) %>%
#   select(tf, id)
#
# chipexo_files = list.files(here('data/chip_exo/raw'), ".tabular$",
#                            recursive = TRUE,
#                            full.names = TRUE)
#
# names(chipexo_files) = str_remove(basename(chipexo_files),
#                                   '_chexmix_allevents.tabular')
#
# read_in_chipexo = function(path){
#
#   meta_line1=read_tsv(path, n_max = 1,skip=1)
#   meta_line2=read_tsv(path, n_max = 1,skip=3)
#
#   read_tsv(path,
#            comment = '#',
#            col_names = c('chr_coord', 'YPD_Sig','YPD_Ctrl','YPD_log2Fold','YPD_log2P','ActiveConds')) %>%
#     separate(chr_coord, into = c('chr','coord'), remove=TRUE, sep=":") %>%
#     mutate(sig_count = meta_line1$TotalSigCount,
#            control_count = meta_line2$CtrlCount,
#            sig_fraction = meta_line1$SignalFraction,
#            sigCtrlScaling = meta_line2$SigCtrlScaling,
#            cond = meta_line1$Name,
#            parent_cond = meta_line2$ParentCond) %>%
#     mutate(YPD_log2P = as.numeric(ifelse(YPD_log2P=='-Infinity', Inf, YPD_log2P)))
#
# }
#
# chipexo_df_list = map(chipexo_files, read_in_chipexo)
#
# chipexo_df = chipexo_df_list %>%
#   bind_rows(.id='chipexo_id') %>%
#   left_join(meta_df) %>%
#   mutate(assay = case_when(
#     assay == "ChIP-exo" ~ 'chipexo',
#     assay == "MNase-ChIP-seq" ~ 'mnase_chip_seq'),
#     tf = toupper(tf)) %>%
#   # CTDSER2 etc is modified RNA Pol II
#   filter(assay == 'chipexo',
#          !tf %in% c('CTDSER2', 'CTDSER5', 'CTDSER7')) %>%
#   left_join(genes_table %>%
#               select(id,gene) %>%
#               dplyr::rename(tf_id=id),
#             by = c('tf' = 'gene')) %>%
#   left_join(chipexo_by_hand_gene_map) %>%
#   mutate(tf_id = ifelse(is.na(tf_id), id, tf_id)) %>%
#   select(-id)

#chipexo_df %>% write_csv("data/chip_exo/chipexo_df.csv.gz")

chipexo_df = read_csv("data/chip_exo/chipexo_df.csv.gz")

# map2(chipexo_df_list, names(chipexo_df_list),
#      ~write_csv(select(.x,chr, coord, YPD_Sig, YPD_Ctrl, YPD_log2Fold, YPD_log2P),
#                 file.path('data/chip_exo/parsed', paste0(.y,'.csv.gz'))))

chipexo_upload_df = chipexo_df %>%
  mutate(file = file.path('data/chip_exo/parsed', paste0(chipexo_id,'.csv.gz'))) %>%
  group_by(chipexo_id, tf_id,replicate, accession,sra_accession,cond,parent_cond,tf, file) %>%
  tally()

chrmap = read_csv("data/genome_files/chrmap.csv.gz")

# map2(chipexo_df_list, names(chipexo_df_list),
#      ~write_csv(select(.x,chr, coord, YPD_Sig, YPD_Ctrl, YPD_log2Fold, YPD_log2P),
#                 file.path('data/chip_exo/parsed', paste0(.y,'.csv.gz'))))

chipexo_upload_df = chipexo_df %>%
  mutate(file = file.path('data/chip_exo/parsed', paste0(chipexo_id,'.csv.gz'))) %>%
  group_by(chipexo_id, tf_id,replicate, accession,sra_accession,cond,parent_cond,tf, file) %>%
  tally() %>%
  left_join(select(genomicfeature_df, id, locus_tag, alias) %>%
              dplyr::rename(tf_locus_tag = locus_tag, tf_alias = alias), by =c('tf_id'='id'))

prepare_http_upload <- function(chipexo_id, tf, tf_locus_tag, cond, parent_cond, file, replicate, ...){

  stopifnot(file.exists(file))

  # read in the file, change coord to start and create a column `end` which
  # is start+1 and save the file to a tmpfile .csv.gz
  tmp <- tempfile(fileext = ".csv.gz")
  read_csv(file) %>%
    mutate(start = coord,
           end = start+1) %>%
    left_join(select(chrmap, numbered, ucsc) %>%
                mutate(numbered = paste0('chr', numbered)),
              by = c('chr' = 'numbered')) %>%
    mutate(chr = ucsc) %>%
    select(-ucsc) %>%
    select(chr,start,end,YPD_Sig,YPD_Ctrl,YPD_log2Fold,YPD_log2P) %>%
    write_csv(tmp)

  # Extract data from the row
  list(
    regulator_locus_tag = tf_locus_tag,
    source_orig_id = chipexo_id,
    condition = cond,
    replicate = replicate,
    file = httr::upload_file(tmp,type="application/gzip"),
    source_name='chipexo_pugh_allevents'
  )
}

url = "http://127.0.0.1:8000/api/binding/"

token = Sys.getenv('local_token')

header <- httr::add_headers(
  "Authorization" = paste("token", token, sep = " ")
)
nrow(chipexo_upload_df)
# note that genomicfeature_df was added to add tf_locus_Tag -- pull that
# from DB. There are at least 3 records which were not added -- check with chipexo_source_id
# from binding_data later
# currently chipexo ids 520 1032 1147
# response_list = map(re_upload, function(i) {
#   message(i)
#   body_data = pmap(chipexo_upload_df[i, ], prepare_http_upload)[[1]]
#   res = POST(url, config = timeout(180), header, body = body_data, encode = "multipart")
#   message(httr::status_code(res))
#   res
# })


binding_df = GET(paste0(url,'export/'), header) %>%
  httr::content() %>%
  filter(source_id == 2)

delete_row <- function(id, ...) {
  message(paste0('deleting id: ', as.character(id)))
  delete_url <- paste0(url, id)
  httr::DELETE(delete_url, header)
}

# map(binding_ids, delete_row)


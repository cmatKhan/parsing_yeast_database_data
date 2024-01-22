library(tidyverse)
library(here)
library(httr)

chrmap = read_csv("data/genome_files/chrmap.csv.gz")

genes_table = read_csv("data/genome_files/genes_table.csv.gz") %>%
  select(-c(uploader,uploadDate,modified,modifiedBy))

cctf = read_csv("data/callingcards/cctf.csv.gz") %>%
  select(-c(uploader,uploadDate,modified,modifiedBy)) %>%
  left_join(genes_table %>% select(id, locus_tag, gene), by = c('tf' = 'id')) %>%
  select(id, tf, locus_tag, gene) %>%
  dplyr::rename(tf_id = tf, tf = id)

ccexperiment = read_csv("data/callingcards/ccexperiment.csv.gz") %>%
  left_join(cctf, by = 'tf') %>%
  select(-c(uploader,uploadDate,modified,modifiedBy))


hops = read_csv("data/callingcards/hops_manifest.csv.gz") %>%
  left_join(ccexperiment, by = c("experiment" = "id")) %>%
  mutate(qbed = file.path("data/callingcards/qbed", basename(qbed))) %>%
  select(-c(uploader,uploadDate,modified,modifiedBy))

# modify_chr = function(chr_format, seqlength, ucsc, qbed, ...){
#   if(chr_format != "ucsc"){
#     data_path = qbed
#
#     df = read_tsv(data_path) %>%
#       left_join(select(chrmap, chr_format, ucsc, seqlength), by =  c('chr' = chr_format)) %>%
#       mutate(end = ifelse(end > seqlength, seqlength, end)) %>%
#       mutate(chr = ucsc) %>%
#       select(-c(ucsc, seqlength))
#
#     df %>% write_tsv(file.path("data/callingcards/ucsc_qbed", basename(data_path)))
#   } else{
#     file.copy(qbed, file.path("data/callingcards/ucsc_qbed", basename(qbed)))
#   }
# }
#
# pmap(hops, modify_chr)

bindingmanualqc = read_csv(here('data/callingcards/manualqc.csv.gz')) %>%
  select(-c(uploader,uploadDate,modified,modifiedBy))

hops = hops %>%
  mutate(qbed = str_replace(qbed,"data/callingcards/qbed", "data/callingcards/ucsc_qbed"),
         chr_format = 'ucsc') %>%
  mutate(qbed = paste0(qbed, '.gz')) %>%
  left_join(bindingmanualqc %>% select(-id) %>% dplyr::rename(qc_notes = note),
            by = c('experiment' = 'experiment'))


prepare_http_upload <- function(id, qbed, locus_tag, batch_replicate,
                                batch, lab, source, chip_better, passing_replicate,
                                data_usable, rank_recall, qc_notes, ...) {

  stopifnot(file.exists(qbed))

  data_source = if(lab == 'mitra'){
    'mitra_cc'
  } else if(lab == 'brent'){
    if(source == 'mitra'){
      'brent_mitra_cc'
    } else{
    "brent_nf_cc"
    }
  }

  best_datatype = if(chip_better == 'no'){
    'pass'
  } else if(chip_better == 'yes'){
    'fail'
  } else{
    chip_better
  }

  passing_replicate = if(passing_replicate == 'yes'){
    'pass'
  } else if(passing_replicate == 'no'){
    'fail'
  } else{
    passing_replicate
  }

  data_usable = if(data_usable == 'yes'){
    'pass'
  } else if(data_usable == 'no'){
    'fail'
  } else{
    data_usable
  }

  # Extract data from the row
  list(
    regulator_locus_tag = locus_tag,
    source_orig_id = id,
    batch = batch,
    replicate = batch_replicate,
    file = httr::upload_file(qbed,type="application/gzip"),
    source_name=data_source,
    rank_recall = rank_recall,
    best_datatype = best_datatype,
    passing_replicate = passing_replicate,
    data_usable = data_usable,
    notes = "source_orig_id is hop id from current cc db",
    qc_notes = qc_notes
  )
}

url = "http://127.0.0.1:8000/api/binding/"

token = Sys.getenv('local_token')

header <- httr::add_headers(
  "Authorization" = paste("token", token, sep = " ")
)
nrow(hops)
# response_list = map(3:nrow(hops), function(i) {
#   message(i)
#   body_data = pmap(hops[i, ], prepare_http_upload)[[1]]
#   res = POST(url, config = timeout(180), header, body = body_data, encode = "multipart")
#   message(httr::status_code(res))
#   Sys.sleep(8)
#   res
# })


binding_df = GET(paste0(url,'export/'), header) %>%
  httr::content() %>%
  filter(source_id == 7 | source_id == 1 | source_id == 8)

delete_row <- function(id, ...) {
  message(paste0('deleting id: ', as.character(id)))
  delete_url <- paste0(url, id)
  httr::DELETE(delete_url, header)
}

responses <- pmap(binding_df, delete_row)

library(tidyverse)

summarize_binding_response_tf = function(tf){

  analysis_expression_df = get_combined_data(
    file.path(url_generator('expression'), 'combined/'),
    list(
      source_time = c('mcisaac_oe,15'),
      regulator_symbol = tf
    ), header) %>%
    left_join(expression_df %>% select(id,name), by = c('record_id' = 'id')) %>%
    dplyr::rename(expression_id = record_id) %>%
    replace_na(list(effect = 1e6, pvalue = 0))

  analysis_binding_df = get_combined_data(
    file.path(url_generator('promotersetsig'), 'combined/'),
    list(
      regulator_symbol = tf
    ), header) %>%
    left_join(promotersetsig_df %>% select(id,name), by = c('record_id' = 'id')) %>%
    dplyr::rename(promotersetsig_id = record_id) %>%
    replace_na(list(effect = 1e6, pvalue = 0))

 x = analysis_binding_df %>%
  dplyr::rename(binding_source = name,
                binding_effect=effect,
                binding_pvalue=pvalue) %>%
  left_join(analysis_expression_df %>%
              dplyr::rename(expression_source = name,
                            expression_effect = effect,
                            expression_pvalue = pvalue), by =c('target_id', 'regulator_symbol'),
            relationship = 'many-to-many') %>%
  mutate(responsive = case_when(
    expression_source == 'mcisaac_oe' & expression_effect > 0 ~ TRUE,
    expression_source == 'kemmeren_tfko' & expression_pvalue < 0.5  ~ TRUE,
    expression_source == 'hu_reimann_tfko' & expression_pvalue < 0.5 ~ TRUE,
    .default = FALSE
  ))

x %>%
  mutate(promotersetsig_id = factor(promotersetsig_id)) %>%
  ungroup() %>%
  filter((!is.na(expression_effect) & !is.na(expression_pvalue))) %>%
  group_by(regulator_symbol, promotersetsig_id, binding_source, expression_source) %>%
  arrange(binding_pvalue) %>%
  slice_head(n=30)

}

# library(tidyverse)
# library(future)
# library(promises)
#
# summarize_binding_response_tf = function(tf) {
#   # Set up future plan for asynchronous execution
#   plan(multisession)
#
#   # Create futures for the asynchronous operations
#   future_analysis_expression = future({
#     get_combined_data(
#       file.path(url_generator('expression'), 'combined/'),
#       list(source_time = c('mcisaac_oe,15'), regulator_symbol = tf), header) %>%
#       left_join(expression_df %>% select(id, name), by = c('record_id' = 'id')) %>%
#       dplyr::rename(expression_id = record_id) %>%
#       replace_na(list(effect = 1e6, pvalue = 0))
#   })
#
#   future_analysis_binding = future({
#     get_combined_data(
#       file.path(url_generator('promotersetsig'), 'combined/'),
#       list(regulator_symbol = tf), header) %>%
#       left_join(promotersetsig_df %>% select(id, name), by = c('record_id' = 'id')) %>%
#       dplyr::rename(promotersetsig_id = record_id) %>%
#       replace_na(list(effect = 1e6, pvalue = 0))
#   })
#
#   # Wait for both futures to resolve
#   analysis_expression_df <- value(future_analysis_expression)
#   analysis_binding_df <- value(future_analysis_binding)
#
#   # Continue with the rest of the function
#   x <- analysis_binding_df %>%
#     dplyr::rename(binding_source = name, binding_effect = effect, binding_pvalue = pvalue) %>%
#     left_join(analysis_expression_df %>%
#                 dplyr::rename(expression_source = name, expression_effect = effect, expression_pvalue = pvalue), by = c('target_id', 'regulator_symbol'),
#               relationship = 'many-to-many') %>%
#     mutate(responsive = case_when(
#       expression_source == 'mcisaac_oe' & expression_effect > 0 ~ TRUE,
#       expression_source == 'kemmeren_tfko' & expression_pvalue < 0.5  ~ TRUE,
#       expression_source == 'hu_reimann_tfko' & expression_pvalue < 0.5 ~ TRUE,
#       .default = FALSE
#     ))
#
#   return(
#     x %>%
#       mutate(promotersetsig_id = factor(promotersetsig_id)) %>%
#       ungroup() %>%
#       filter((!is.na(expression_effect) & !is.na(expression_pvalue))) %>%
#       group_by(regulator_symbol, promotersetsig_id, binding_source, expression_source) %>%
#       arrange(binding_pvalue) %>%
#       slice_head(n = 30)
#   )
# }
#
#

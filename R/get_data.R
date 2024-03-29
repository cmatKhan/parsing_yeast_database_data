library(httr)
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)
library(future)
library(furrr)

get_data <- function(url, token, num_records = 10) {
  # Set up future for parallel processing
  plan(multisession)

  header <- httr::add_headers(
    "Authorization" = paste("token", token, sep = " "),
    "Content-Type" = "application/json"
  )

  # Initial request to get total record count
  response <- GET(url, query = list(limit = 1), header)
  if(status_code(response) != 200) {
    stop("Failed to fetch data: ", status_code(response))
  }
  total_records <- content(response, "parsed")$count

  # Calculate the number of chunks
  num_chunks <- ceiling(total_records / num_records)
  message(paste0('fetching: ', as.character(num_chunks)))

  # Function to fetch a chunk of data
  fetch_chunk <- function(chunk_id, num_records) {
    message(paste0('fetching chunk: ', as.character(chunk_id)))
    offset <- (chunk_id - 1) * num_records
    params <- list(limit = num_records, offset = offset)
    response <- GET(url, query = params)
    if(status_code(response) != 200) {
      return(tibble())  # Return an empty tibble in case of error
    }
    data <- content(response, "parsed")
    gene_data <- data$results
    return(map_dfr(gene_data, as_tibble))
  }

  # Fetch all chunks in parallel
  results <- future_map_dfr(1:num_chunks, ~fetch_chunk(.x, num_records))

  return(results)
}

# Example usage
# url <- "http://3.143.144.93/api/v1/genes/"
# token = 'lakjsfdlkasjdf'
# result <- get_data(url, token, num_records = 500)


recursive_get <- function(url, header, res_vector = list()) {
  res <- httr::GET(url, header)

  # Append the current response to the result vector
  res_vector <- c(res_vector, list(res))

  # Extract 'next' URL from the response
  next_url <- httr::content(res)$`next`

  # Recursive call if 'next' URL exists
  if (!is.null(next_url)) {
    res_vector <- recursive_get(next_url, header, res_vector)
  }

  return(res_vector)
}

# res = recursive_get(url,header)
# in example of binding endpoint, `file` was sometimes NULL
# for each result object, process prior to transforming to
# dataframe

process_binding_results = function(res) {
  map(httr::content(res)$`results`, ~{
    if(is.null(.x$file)) {
      .x$file <- ''
    }
    .x
  }) %>%
    map_dfr(as_tibble)
}

# Example usage
# bind_rows(map(res, process_binding_results))

# to delete all rows
# delete_row = function(id,...){
#   delete_url = paste0(url,id,'/')
#   message(delete_url)
#   httr::DELETE(delete_url,header)
# }

library(httr)
library(dplyr)
library(tidyr)
library(purrr)
library(jsonlite)
library(future)
library(furrr)

get_data <- function(url, num_records = 10) {
  # Set up future for parallel processing
  plan(multisession)

  # Initial request to get total record count
  response <- GET(url, query = list(limit = 1))
  if(status_code(response) != 200) {
    stop("Failed to fetch data: ", status_code(response))
  }
  total_records <- content(response, "parsed")$count

  # Calculate the number of chunks
  num_chunks <- ceiling(total_records / num_records)

  # Function to fetch a chunk of data
  fetch_chunk <- function(chunk_id, num_records) {
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
# result <- get_data(url, num_records = 500)

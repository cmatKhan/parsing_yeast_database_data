library(httr)
library(jsonlite)

get_data_file <- function(url, token, params_list, directory_path) {
  # Make the API request with the provided parameters
  header <- httr::add_headers(
    "Authorization" = paste("token", token, sep = " "),
    "Content-Type" = "application/json"
  )
  response <- GET(url, query = params_list,header)
  if (status_code(response) != 200) {
    stop("Failed to fetch data: ", status_code(response))
  }

  # Parse the response
  data <- content(response, "parsed")

  # Extract the file URL from the response
  file_url <- data$results[[1]]$file
  if (is.null(file_url)) {
    stop("No file URL found in the response")
  }

  # Extract the basename of the file URL
  file_name <- basename(file_url)

  # Construct the file path for saving
  file_path <- file.path(directory_path, file_name)

  # Download and save the file
  download.file(file_url, file_path, mode = "wb")

  cat("File downloaded and saved to", file_path, "\n")
  return(invisible())
}

# Example usage
# url <- "http://example.com/api/data"
# params_list <- list(experiment_id=32, promoter_source='yiming', background_source='adh1', hops_source='nf-core/callingcards:dev')
# get_data_file(url, params_list, directory_path = "/path/to/directory")


# Example usage
# url <- "http://3.143.144.93/api/v1/genes/"
# params_list=list(source='nf-core/callingcards:dev', tf_locus_tag='HAP4')
# x = tempdir()
# get_data_file(url,params_list,x)

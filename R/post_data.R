library(dplyr)
library(jsonlite)
library(httr)
library(purrr)


# Function to format error message
format_error_message <- function(response) {
  content <- httr::content(response)
  error_message <- "Error: "
  for (field in names(content)) {
    errors <- content[[field]]
    for (error in errors) {
      error_message <- paste0(error_message, field, " - ", error, "; ")
    }
  }

  return(error_message)
}


# Define the function with token authentication
post_data <- function(tibble_data, url, token) {
  # Prepare the header for token authentication
  header <- httr::add_headers(
    "Authorization" = paste("token", token, sep = " "),
    "Content-Type" = "application/json"
  )

  # Function to POST a single row
  post_row <- function(i) {
    row_data = tibble_data[i,]
    # Convert the row to JSON
    json_data <- toJSON(as.list(row_data), auto_unbox = TRUE)

    # POST the JSON to the database with the authentication header
    response <- POST(url, body = json_data, encode = "json", config = header)

    # Optional: Check response status
    if (status_code(response) != 201) {
      warning(paste("Failed to POST data at row", i, ":", content(response, "text", encoding = "UTF-8")))
    }

    response
  }

  # Apply the function to each row
  map(1:nrow(tibble_data), post_row)
}

# Example usage
# Assume `my_tibble` is your tibble, `database_url` is the URL of your database, and `my_token` is your authentication token
# post_data(my_tibble, database_url, my_token)

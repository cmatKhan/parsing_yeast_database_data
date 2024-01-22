library(tidyverse)
library(httr)

get_combined_data = function(url, query, header){
  res = GET(url, query=query, header)

  if (httr::status_code(res) != 200) {
    message("error: ", httr::status_code(res))
    message(as.character(query))
  } else{

    temp_gz <- tempfile(fileext = ".gz")

    writeBin(res %>% httr::content(), temp_gz)

    df = read_csv(temp_gz)

    file.remove(temp_gz)

    return(df)
  }


}

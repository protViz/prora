#' map_ids_uniprot
#'
#' Access to uniprot web API. For more detail see https://www.uniprot.org/help/uploadlists
#'
#' @param data data
#' @param ID_col column name containing id's
#' @param from id type
#' @param to id type
#' @param format in which format to return the mapping
#' @importFrom httr POST content
#' @importFrom rlang :=
#' @export
#'
#' @examples
#' library(tidyverse)
#' fc_estimates <- prora::exampleContrastData
#'
#' filtered_dd <- get_UniprotID_from_fasta_header(fc_estimates, idcolumn = "protein_Id")
#'
#' map_ids_uniprot( filtered_dd )
#'
#'
map_ids_uniprot <- function(data,
                            ID_col = "UniprotID",
                            from =  "ACC+ID",
                            to = "P_ENTREZGENEID",
                            format = "tab"){

  ids_to_map <- unique(na.omit(data[[ID_col]]))

  url = "https://www.uniprot.org/uploadlists/"
  params = list(
    from = from,
    to = to,
    format = format,
    query = paste(ids_to_map, collapse = " ")
  )

  r <- httr::POST(url, body = params, encode = "form")
  mapping <- readr::read_tsv(httr::content(r))
  mapping <- mapping %>% dplyr::rename( !!ID_col := "From", !!to := "To" )
  res <- dplyr::right_join(mapping, data, by = ID_col)

  return(res)
}

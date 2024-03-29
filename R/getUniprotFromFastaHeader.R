#' Translate Fasta header to UniprotSwissprot ID
#'
#' @export
#'
#' @param df \code{data.frame} with FASTA headers in first column
#' @param idcolumn Column name containing the IDs
#'
#' @return Returns the whole \code{data.frame} with a column called \code{UniprotID} containing UniprotSwissprot IDs instead of Fasta headers. This format is easier to use for subsequent ID mappings using the \code{Annotationdbi} package.
#'
#' @examples
#' data("exampleContrastData", package = "prora")
#' get_UniprotID_from_fasta_header(exampleContrastData)
get_UniprotID_from_fasta_header <- function(df, idcolumn = "protein_Id") {
  map <- df %>% dplyr::select(idcolumn) %>% distinct() %>%
    dplyr::filter(grepl(pattern = "^sp|^tr", !!sym(idcolumn))) %>%
    tidyr::separate(col = !!sym(idcolumn),
                    sep = "_",
                    into = c("begin", "end"),
                    remove = FALSE) %>%
    tidyr::separate(
      col = !!sym("begin"),
      sep = "\\|",
      into = c("prefix", "UniprotID", "Symbol")
    ) %>%
    dplyr::select(-!!sym("prefix"), -!!sym("Symbol"), -!!sym("end"))
  res <- dplyr::right_join(map, df , by= idcolumn)
  return(res)
}

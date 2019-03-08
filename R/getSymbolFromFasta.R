#' Translate Fasta header to gene Symbol
#'
#' @export
#' @importFrom dplyr filter select
#' @importFrom tidyr separate
#'
#' @param df \code{data.frame} with FASTA headers in first column
#'
#' @return Returns the whole \code{data.frame} with a column called \code{Symbol} containing gene symbols instead of Fasta headers. This format is easier to use for subsequent ID mappings using the \code{Annotationdbi} package.
#'
#' @examples
#' data("exampleContrastData", package = "fgczgseaora")
#' getSymbolFromSwissprotID(exampleContrastData)
getSymbolFromSwissprotID <- function(df) {
  colnames(df)[1] <- "protein_Id"
  df %>%
    dplyr::filter(grepl(pattern = "sp", df$protein_Id)) %>%
    tidyr::separate(col = protein_Id,
             sep = "_",
             into = c("begin", "end")) %>%
    tidyr::separate(
      col = begin,
      sep = "\\|",
      into = c("prefix", "uniprotid", "Symbol")
    ) %>%
    dplyr::select(-prefix, -uniprotid, -end) %>%
    return()
}

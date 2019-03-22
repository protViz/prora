#' Translate Fasta header to UniprotSwissprot ID
#'
#' @export getUniprotFromFastaHeader
#' @importFrom dplyr filter select
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#'
#' @param df \code{data.frame} with FASTA headers in first column
#'
#' @return Returns the whole \code{data.frame} with a column called \code{UniprotID} containing UniprotSwissprot IDs instead of Fasta headers. This format is easier to use for subsequent ID mappings using the \code{Annotationdbi} package.
#'
#' @examples
#' data("exampleContrastData", package = "fgczgseaora")
#' getUniprotFromFastaHeader(exampleContrastData)
getUniprotFromFastaHeader <- function(df) {
  colnames(df)[1] <- "protein_Id"
  df %>%
    dplyr::filter(grepl(pattern = "sp", df$protein_Id)) %>%
    tidyr::separate(col = protein_Id,
             sep = "_",
             into = c("begin", "end")) %>%
    tidyr::separate(
      col = begin,
      sep = "\\|",
      into = c("prefix", "UniprotID", "Symbol")
    ) %>%
    dplyr::select(-prefix, -Symbol, -end) %>%
    return()
}

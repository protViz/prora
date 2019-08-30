#' Translate Fasta header to UniprotSwissprot ID
#'
#' @export getUniprotFromFastaHeader
#' @importFrom dplyr filter select
#' @importFrom tidyr separate
#' @importFrom magrittr %>%
#' @importFrom rlang sym
#'
#' @param df \code{data.frame} with FASTA headers in first column
#' @param idcolumn Column name containing the IDs
#'
#' @return Returns the whole \code{data.frame} with a column called \code{UniprotID} containing UniprotSwissprot IDs instead of Fasta headers. This format is easier to use for subsequent ID mappings using the \code{Annotationdbi} package.
#'
#' @examples
#' data("exampleContrastData", package = "fgczgseaora")
#' getUniprotFromFastaHeader(exampleContrastData)
getUniprotFromFastaHeader <- function(df, idcolumn = "protein_Id") {
  df %>%
    dplyr::filter(grepl(pattern = "sp", !!sym(idcolumn))) %>%
    tidyr::separate(col = !!sym(idcolumn),
                    sep = "_",
                    into = c("begin", "end")) %>%
    tidyr::separate(
      col = !!sym("begin"),
      sep = "\\|",
      into = c("prefix", "UniprotID", "Symbol")
    ) %>%
    dplyr::select(-!!sym("prefix"), -!!sym("Symbol"), -!!sym("end")) %>%
    return()
}

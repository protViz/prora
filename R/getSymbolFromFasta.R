#' Translate Fasta header to gene Symbol
#'
#' @export getSymbolFromFasta getSymbolFromFasta
#' @import tidyverse sigora org.Hs.eg.db GO.db reactome.db dplyr AnnotationDbi
#'
#' @param df \code{data.frame} with column called \code{protein_Id}
#'
#' @return Returns the whole \code{data.frame} with a column called \code{Symbol} containing gene symbols instead of Fasta headers. This format is easier to use for subsequent ID mappings using the \code{Annotationdbi} package.
#'
#' @examples
#' data("exampleContrastData", package = "fgczgseaora")
#' getSymbolFromFasta(exampleContrastData)

getSymbolFromFasta <- function(df) {
  df %>%
    dplyr::filter(grepl(pattern = "sp", df$protein_Id)) %>%
    separate(col = protein_Id,
             sep = "_",
             into = c("begin", "end")) %>%
    separate(
      col = begin,
      sep = "\\|",
      into = c("prefix", "uniprotid", "Symbol")
    ) %>%
    dplyr::select(-prefix, -uniprotid, -end) %>%
    return()
}

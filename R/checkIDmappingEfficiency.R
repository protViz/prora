#' Checks mapping efficiency for IDs on the most commonly used ID types
#' @param IDs character vector of IDs
#' @param keytype corresponding keytype
#' @importFrom AnnotationDbi mapIds
#' @importFrom S4Vectors na.omit
#' @importFrom dplyr distinct
#' @importFrom tibble enframe
#' @return Returns a data.frame containing information of mapping efficiency on GO, ENTREZ and KEGG IDs
#'
#' @examples
#' data("exampleContrastData", package = "fgczgseaora")
#' df <- getUniprotFromFastaHeader(exampleContrastData)
#' checkIDmappingEfficiency(df$UniprotID, keytype = "UNIPROT")
#'
#' @export checkIDmappingEfficiency
#'
checkIDmappingEfficiency <- function(IDs, keytype) {
  uniqueInputIDs <- unique(IDs)
  number_unique_InputIDs <- length(uniqueInputIDs)

  mapping <- function(x = org.Hs.eg.db, column) {
    mpd <- mapIds(
      x = x,
      keys = uniqueInputIDs,
      keytype = keytype,
      column = column,
      multiVals = "CharacterList"
    ) %>% unlist %>%
      enframe(name = "InputIDs", value = "mappedIDs") %>%
      na.omit %>% distinct(InputIDs) %>% nrow %>% return()
  }

  ENTREZmapped <- mapping(column = "ENTREZID")
  KEGGmapped <- mapping(column = "PATH")
  GOmapped <- mapping(column = "GO")

  out <- data.frame(InputIDs = c(number_unique_InputIDs, 100),
       ENTREZIDs = c(ENTREZmapped, ENTREZmapped/number_unique_InputIDs*100),
       KEGGIDs = c(KEGGmapped, KEGGmapped/number_unique_InputIDs*100),
       GOIDs = c(GOmapped, GOmapped/number_unique_InputIDs*100)
       )
  row.names(out) <- c("Number", "Percent Mapped")
  return(out)
}

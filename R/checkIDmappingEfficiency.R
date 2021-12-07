#' Checks mapping efficiency for IDs on the most commonly used ID types
#' @param IDs character vector of IDs
#' @param keytype corresponding keytype
#' @param db database object, default org.Hs.eg.db
#' @return Returns a data.frame containing information of mapping efficiency on GO, ENTREZ and KEGG IDs
#'
#'
#' @export checkIDmappingEfficiency
#'
#' @examples
#' library(org.Hs.eg.db)
#' data("exampleContrastData", package = "prora")
#' df <- get_UniprotID_from_fasta_header(exampleContrastData)
#' checkIDmappingEfficiency(df$UniprotID, keytype = "UNIPROT")
#'
checkIDmappingEfficiency <- function(IDs, keytype, db = org.Hs.eg.db) {
  uniqueInputIDs <- unique(IDs)
  number_unique_InputIDs <- length(uniqueInputIDs)

  mapping <- function(column) {
    mpd <- AnnotationDbi::mapIds(
      x = db,
      keys = uniqueInputIDs,
      keytype = keytype,
      column = column,
      multiVals = "CharacterList"
    ) %>% unlist %>%
      enframe(name = "InputIDs", value = "mappedIDs") %>%
      na.omit %>% distinct(!!sym("InputIDs")) %>% nrow %>% return()
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

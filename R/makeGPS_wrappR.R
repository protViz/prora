#' makeGPS Wrapper
#'
#' Function for generating background GPS repository for sigora and ora
#'
#' @export makeGPS_wrappR makeGPS_wrappR
#'
#' @import tidyverse sigora org.Hs.eg.db GO.db reactome.db dplyr AnnotationDbi
#'
#'
#' @param ids Character vector of gene symbols (experimental background)
#' @param target Character, database to be used to generate GPS repository ("KEGG", "GO", "reactome")
#' @return A GPS repository based on the self specified background in \code{ids}, see \code{sigora::makeGPS()} for more information.
#'
#' @seealso \code{sigora::makeGPS()}
#'
#' @examples
#' data("exampleSymbols", package = "fgczgseaora")
#' myGPSrepo <- makeGPS_wrappR(ids = exampleSymbols, target = "KEGG")
#'

makeGPS_wrappR <- function(ids, target = "KEGG") {
  if (target == "KEGG") {
    gp_db <- org.Hs.eg.db
    target_column <- "PATH"
    pn_table <- sigora::kegH$pathwaydescriptions
    colnames(pn_table) <- c("pathwayID", "pathwayName")
    pn_table$pathwayID <-
      substr(pn_table$pathwayID, start = 4, stop = 8)
    gp_table <- AnnotationDbi::mapIds(
      gp_db,
      keys = ids,
      keytype = "SYMBOL",
      column = target_column,
      multiVals = "CharacterList"
    ) %>% unlist %>% data.frame(Symbol = names(.), pathwayID = .)
  } else {
    if (target == "reactome") {
      gp_db <- reactome.db
      pn_db <- reactome.db
      ids <- AnnotationDbi::mapIds(
        # First map to ENTREZID
        org.Hs.eg.db,
        keys = ids,
        keytype = "SYMBOL",
        column = "ENTREZID",
        multiVals = "CharacterList"
      ) %>% unlist %>% na.omit
      target_column <- "PATHID"
      k0 = "ENTREZID"
      k1 = "PATHID"
      k2 = "PATHNAME"
    } else
      if (target == "GO") {
        gp_db <- org.Hs.eg.db
        pn_db <- GO.db
        target_column = "GO"
        k0 = "SYMBOL"
        k1 = "GOID"
        k2 = "TERM"
      } else
        return(message("Specify a valid target database (KEGG, reactome, GO)"))

    gp_table <- mapIds(
      gp_db,
      keys = ids,
      keytype = k0,
      column = target_column,
      multiVals = "CharacterList"
    ) %>% unlist %>% data.frame(Symbol = names(.), pathwayID = .)

    pn_table <- mapIds(
      pn_db,
      keys = as.character(gp_table$pathwayID),
      multiVals = "CharacterList",
      keytype = k1,
      column = k2
    ) %>% unlist %>% data.frame(pathwayID = names(.), pathwayName = .)
  }

  mkGPStable <- inner_join(pn_table, gp_table) %>% distinct

  colnames(mkGPStable) <- c("pathwayId", "pathwayName", "gene")

  out <- makeGPS(mkGPStable)

  return(out)
}

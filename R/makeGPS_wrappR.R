#' makeGPS Wrapper
#' Function for generating background GPS repository for sigora and ora
#'
#' @export makeGPS_wrappR makeGPS_wrappR
#'
#' @importFrom org.Hs.eg.db org.Hs.eg.db
#' @importFrom GO.db GO.db
#' @importFrom reactome.db reactome.db
#' @importFrom dplyr inner_join
#' @importFrom sigora makeGPS
#' @importFrom S4Vectors na.omit
#' @importFrom AnnotationDbi mapIds
#' @importFrom magrittr %>%
#' @param ids Character vector of gene symbols (experimental background)
#' @param target Character, database to be used to generate GPS repository ("KEGG", "GO", "reactome")
#' @param dev Logical, if \code{TRUE} the function returns the \code{data.frame} \code{makeGPS} needs as input
#' @return A GPS repository based on the self specified background in \code{ids}, see \code{\link[sigora]{makeGPS}} for more information.
#'
#' @seealso \code{\link[sigora]{makeGPS}}
#'
#' @examples
#'
#' data("exampleUniprotIDs", package = "fgczgseaora")
#' myGPSrepo <- makeGPS_wrappR(ids = exampleUniprotIDs)
#' myGPSrepo <- makeGPS_wrappR(ids = exampleUniprotIDs, target = "GO")
#' myGPSrepo <- makeGPS_wrappR(ids = exampleUniprotIDs, target = "react")
makeGPS_wrappR <- function(ids, target = c("KEGG","GO","reactome"), dev = FALSE) {
  target <- match.arg(target)
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
      keytype = "UNIPROT",
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
        keytype = "UNIPROT",
        column = "ENTREZID",
        multiVals = "CharacterList"
      ) %>% unlist %>% S4Vectors::na.omit()
      target_column <- "PATHID"
      k0 = "ENTREZID"
      k1 = "PATHID"
      k2 = "PATHNAME"
    } else
      if (target == "GO") {
        gp_db <- org.Hs.eg.db
        pn_db <- GO.db
        target_column = "GO"
        k0 = "UNIPROT"
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

  mkGPStable <- dplyr::inner_join(pn_table, gp_table) %>% dplyr::distinct()

  colnames(mkGPStable) <- c("pathwayId", "pathwayName", "gene")

  if(dev) return(mkGPStable)

  out <- makeGPS(mkGPStable)

  return(out)
}

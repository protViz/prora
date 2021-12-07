#' makeGPS Wrapper
#' Function for generating background GPS repository for sigora and ora
#'
#' @export makeGPS_wrappR makeGPS_wrappR
#'
#' @param ids Character vector of gene symbols (experimental background)
#' @param target Character, database to be used to generate GPS repository ("KEGG", "GO", "reactome")
#' @param keytype default UNIPROT
#' @return A GPS repository based on the self specified background in \code{ids}, see \code{\link[sigora]{makeGPS}} for more information.
#'
#' @seealso \code{\link[sigora]{makeGPS}}
#'
#' @examples
#'
#' data("exampleUniprotIDs", package = "prora")
#' myGPSrepo <- makeGPS_wrappR(ids = exampleUniprotIDs)
#' \dontrun{
#' myGPSrepo <- makeGPS_wrappR(ids = exampleUniprotIDs, target = "GO")
#' myGPSrepo <- makeGPS_wrappR(ids = exampleUniprotIDs, target = "react")
#' }
makeGPS_wrappR <- function(ids, keytype = "UNIPROT",
                           target = c("KEGG", "GO", "reactome")) {
  target <- match.arg(target)

  if (target == "KEGG") {
    tabs <- .get_tables_KEGG(ids = ids, keytype = "UNIPROT")
  } else if (target %in% c("GO", "reactome")) {
    tabs <- .get_tables_ractome_GO(ids = ids, keytype = "UNIPROT",  target = target)
  } else {
    stop("Specify a valid target database (one of KEGG, GO, reactome)")
  }

  mkGPStable <-
    dplyr::inner_join(tabs$pn_tab, tabs$gp_tab) %>%
    dplyr::distinct()

  colnames(mkGPStable) <- c("pathwayId", "pathwayName", "gene")

  if (all(dim(mkGPStable) > 0)) {
    return(list(gps = sigora::makeGPS(as.data.frame(mkGPStable)),gpsTable = mkGPStable))
  } else {
    stop("No pathway, gene, description mappings. Cannot produce GPS repository.")
  }
}

.get_tables_ractome_GO <- function(ids, keytype = "UNIPROT", target = c("reactome",  "GO")) {
  target <- match.arg(target)
  if (target == "reactome") {
    gp_db <- reactome.db::reactome.db
    pn_db <- reactome.db::reactome.db

    ids <- AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = ids,
      keytype = keytype,
      column = "ENTREZID",
      multiVals = "CharacterList"
    ) %>% unlist %>% S4Vectors::na.omit()

    target_column <- "PATHID"
    k0 = "ENTREZID"
    k1 = "PATHID"
    k2 = "PATHNAME"
  } else if (target == "GO") {
    library(org.Hs.eg.db)
    gp_db <- org.Hs.eg.db
    library(GO.db)
    pn_db <- GO.db

    target_column = "GO"
    k0 = "UNIPROT"
    k1 = "GOID"
    k2 = "TERM"
  }

  gp_table <- mapIds(
    gp_db,
    keys = ids,
    keytype = k0,
    column = target_column,
    multiVals = "CharacterList"
  ) %>% unlist %>% enframe(name = "Symbol", value = "pathwayID")

  pn_table <- mapIds(
    pn_db,
    keys = as.character(gp_table$pathwayID),
    multiVals = "CharacterList",
    keytype = k1,
    column = k2
  ) %>% unlist %>% enframe(name = "pathwayID", value = "pathwayName")


  out <- list(gp_tab = gp_table, pn_tab = pn_table)
  return(out)
}

.get_tables_KEGG <- function(ids, keytype = "UNIPROT") {
  gp_db <- org.Hs.eg.db
  target_column <- "PATH"

  pn_table <- sigora::load_data('kegH')$pathwaydescriptions %>%
    dplyr::rename(pathwayID = !!sym("pwys"), pathwayName = !!sym("nms")) %>%
    mutate(pathwayID = substr(!!sym("pathwayID"), start = 4, stop = 8))

  gp_table <- AnnotationDbi::mapIds(
    gp_db,
    keys = ids,
    keytype = keytype,
    column = target_column,
    multiVals = "CharacterList"
  ) %>% unlist %>% enframe(name = "Symbol", value = "pathwayID")

  out <- list(pn_tab = pn_table, gp_tab = gp_table)
  return(out)
}

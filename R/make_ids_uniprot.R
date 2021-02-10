#' map_ids_uniprot
#'
#' Access to uniprot web API. For more detail see https://www.uniprot.org/help/uploadlists
#'
#' @param data data
#' @param ID_col column name containing id's
#' @param from id type
#' @param to id type
#' @param format in which format to return the mapping
#' @importFrom httr POST content
#' @importFrom rlang :=
#' @export
#'
#' @examples
#' library(tidyverse)
#' fc_estimates <- prora::exampleContrastData
#'
#' filtered_dd <- get_UniprotID_from_fasta_header(fc_estimates, idcolumn = "protein_Id")
#' debug(map_ids_uniprot)
#' map_ids_uniprot( filtered_dd )
#'
#'
map_ids_uniprot <- function(data,
                            ID_col = "UniprotID",
                            from =  "ACC+ID",
                            to = "P_ENTREZGENEID",
                            format = "tab"){

  ids_to_map <- unique(na.omit(data[[ID_col]]))

  url = "https://www.uniprot.org/uploadlists/"
  params = list(
    from = from,
    to = to,
    format = format,
    query = paste(ids_to_map, collapse = " ")
  )

  r <- httr::POST(url, body = params, encode = "form")
  bb <- httr::content(r)
  class(bb)
  mapping <- readr::read_tsv(bb)
  mapping <- mapping %>% dplyr::rename( !!ID_col := "From", !!to := "To" )
  res <- dplyr::right_join(mapping, data, by = ID_col)

  return( res )
}

.get_species_mapping_AnnotationHub <- function(species = c("Homo sapiens", "Mus musculus")) {
  species <- match.arg(species)

  ah <- AnnotationHub::AnnotationHub()
  orgdb <- AnnotationHub::query(ah, c("OrgDb", "maintainer@bioconductor.org"))
  if ( !species %in% orgdb$species) {
    stop(paste0("species not found : ", species, "\n"))
  }
  specODB <- orgdb[[which(species == orgdb$species)]]
  egid <- AnnotationDbi::keys(specODB, "ENTREZID")
  res <- AnnotationDbi::select(specODB, egid, c("SYMBOL", "GENENAME", "UNIPROT"), "ENTREZID")
  return(res)

}

#' map id Annotation Hub
#'
#' @export
#' @examples
#'
#' library(tidyverse)
#' fc_estimates <- prora::exampleContrastData
#' #fc_estimates %>% filter(!(grepl("^REV_", protein_Id) | grepl("^zz", protein_Id) | grepl("^CON__", protein_Id)))
#' fc_estimates <- fc_estimates %>% filter(!(grepl("^REV_|^zz|^CON__", protein_Id)))
#' filtered_dd <- get_UniprotID_from_fasta_header(fc_estimates, idcolumn = "protein_Id")
#' dim(filtered_dd)
#' tmp <- map_ids_annotationHub(filtered_dd, species = "Homo sapiens")
#' sum(is.na(tmp$P_ENTREZGENEID ))
#'
map_ids_annotationHub <- function(x, ID_col = "UniprotID", species =  c("Homo sapiens", "Mus musculus")){
  species <- match.arg(species)
  res <- .get_species_mapping_AnnotationHub(species)
  res <- res %>% dplyr::select(!!ID_col := UNIPROT, P_ENTREZGENEID = ENTREZID  )
  res <- na.omit(res)
  res <- right_join(res,  x,  by = ID_col)
  return(res)
}

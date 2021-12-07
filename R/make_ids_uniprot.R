#' map_ids_uniprot
#'
#' Access to uniprot web API. For more detail see https://www.uniprot.org/help/uploadlists
#' @family id_mapping
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
#' # debug(map_ids_uniprot)
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

  if (Sys.info()["sysname"] == "Windows") {
    r <- httr::POST(url, body = params, encode = "form")
  } else {
    # This code was added to handle SSL issues
    # https://github.com/r-lib/httr/issues/669
    httr_config <- httr::config(ssl_cipher_list = "DEFAULT@SECLEVEL=1")
    r <- httr::with_config(config = httr_config, httr::POST(url, body = params, encode = "form"))
  }

  #r <- httr::POST(url, body = params, encode = "form")
  bb <- httr::content(r)
  class(bb)
  mapping <- readr::read_tsv(bb)
  mapping <- mapping %>% dplyr::rename( !!ID_col := "From", !!to := "To" )
  mapping[[to]] <- as.character(mapping[[to]])
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
#' @family id_mapping
#' @param x data frame
#' @param ID_col column with ID's
#' @param species default "Homo sapiens"
#' @export
#' @examples
#'
#'
#' library(tidyverse)
#' fc_estimates <- prora::exampleContrastData
#' \dontrun{
#' fc_estimates %>%
#'  filter(!(grepl("^REV_", protein_Id) | grepl("^zz", protein_Id) | grepl("^CON__", protein_Id)))
#' }
#' fc_estimates <- fc_estimates %>% filter(!(grepl("^REV_|^zz|^CON__", protein_Id)))
#' filtered_dd <- get_UniprotID_from_fasta_header(fc_estimates, idcolumn = "protein_Id")
#' dim(filtered_dd)
#' tmp <- map_ids_annotationHub(filtered_dd, species = "Homo sapiens")
#' sum(is.na(tmp$P_ENTREZGENEID ))
#'
map_ids_annotationHub <- function(x, ID_col = "UniprotID", species =  c("Homo sapiens", "Mus musculus")){
  species <- match.arg(species)
  res <- .get_species_mapping_AnnotationHub(species)
  res <- res %>% dplyr::select(!!ID_col := .data$UNIPROT, P_ENTREZGENEID = .data$ENTREZID  )
  res <- na.omit(res)
  res <- right_join(res,  x,  by = ID_col)
  return(res)
}



#' map ids two ways
#' @family id_mapping
#' @export
#' @param clusterAssignment data.frame
#' @param ID_col column with IDS
#' @param species default "Homo sapiens"
map_ids_2ways <- function(clusterAssignment,
                          ID_col = "UniprotID",
                          species =  c("Homo sapiens", "Mus musculus"))
{
  # map to entriz id's
  .ehandler = function(e) {
    warning("WARN :", e)
    # return string here
    as.character(e)
  }
  #undebug(prora::map_ids_annotationHub)
  result <- tryCatch(
    prora::map_ids_uniprot(clusterAssignment, ID_col = ID_col),
    error = .ehandler)

  id.mapping.service <- "UniProt"
  if (is.character(result)) {
    result <- tryCatch(
      prora::map_ids_annotationHub(clusterAssignment, ID_col = ID_col, species = species),
      error = .ehandler)
    id.mapping.service <- "AnnotationHub"
  }

  return( list(clusterAssignment = result, mapping.service = id.mapping.service) )

}



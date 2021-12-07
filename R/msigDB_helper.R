#' getMsigdbGenesets
#' @family fgsea
#' @param msigCollection data.frema with columns gs_cat, gs_subcat
#' @param species e.g. "Homo sapiens"
#' @export
#' @examples
#' library(tidyverse)
#' library(msigdbr)
#' msigdbr::msigdbr_species()
#' species <- "Homo sapiens"
#' species <- "Mus musculus"
#'
#' hallmark <- msigdbr_collections() %>% filter(.data$gs_cat == "H")
#'
#' #hallmark$gs_subcat <- "HALLMARK"
#' C5 <- bind_rows( {msigdbr_collections() %>%
#'  filter(.data$gs_cat == "C5") %>%
#'  filter(grepl("^GO:", .data$gs_subcat))},
#'  hallmark,
#'  {msigdbr_collections() %>% filter(.data$gs_subcat == "CP:KEGG")} )
#'
#' C5
#' fgseaGSlist <- prora::getMsigdbGenesets(C5, species)
#'
getMsigdbGenesets <- function(msigCollection, species){
  genesetsC5 <- vector(mode = "list", length = nrow(msigCollection) )
  for (i in 1:nrow(msigCollection)) {
    genesetsC5[[i]] <- msigdbr::msigdbr(species = species,
                               category = msigCollection$gs_cat[i],
                               subcategory = msigCollection$gs_subcat[i])
  }
  names(genesetsC5) <- make.names(msigCollection$gs_subcat)
  fgseaGSlist <- lapply(genesetsC5 , fgsea_msigdb)
  return(fgseaGSlist)
}

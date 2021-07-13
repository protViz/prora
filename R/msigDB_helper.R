#' getMsigdbGenesets
#' @export
#' @examples
#' #hallmark$gs_subcat <- "HALLMARK"
#' #library(msigdbr)
#'
#' #C5 <- dplyr::bind_rows( msigdbr::msigdbr_collections() %>% filter(gs_cat == "C5") %>% filter(grepl("^GO:", gs_subcat)),
#' #                 hallmark,
#' #                 msigdbr_collections() %>% filter(gs_subcat == "CP:KEGG") )
#' #species <- "Mus musculus"
getMsigdbGenesets <- function(msigCollection, species){
  genesetsC5 <- vector(mode = "list", length = nrow(msigCollection) )
  for (i in 1:nrow(msigCollection)) {
    genesetsC5[[i]] <- msigdbr(species = species,
                               category = msigCollection$gs_cat[i],
                               subcategory = msigCollection$gs_subcat[i])
  }
  names(genesetsC5) <- make.names(msigCollection$gs_subcat)
  fgseaGSlist <- lapply(genesetsC5 , fgsea_msigdb)
  return(fgseaGSlist)
}

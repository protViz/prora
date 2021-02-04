#' get fgsea compatible rank list from data.frame
#' @export
fgsea_rank <- function(df,
                       ids = "X1",
                       score = "X2"){
  ranks <- df[[score]]
  names(ranks) <- df[[ids]]
  return(ranks)
}


#' convert data frame into list of ranks
#' @export
fgsea_rank_contrasts <-
  function(df,
           ids = "X1",
           score = "X2",
           contrast = "contrast"){
    df <- select(df, dplyr::all_of(c(contrast, ids, score)))
    df <- na.omit(df)
    ldf <-  dplyr::group_by(df, !!sym(contrast)) %>% nest()

    allrnk <- vector(mode = "list", length = nrow(ldf))
    for (i in 1:nrow(ldf)) {
      ranks <- fgsea_rank(ldf$data[[i]], ids = ids, score = score)
      allrnk[[i]] <- ranks
    }
    names(allrnk) <- ldf[[1]]
    return(allrnk)
  }


#' converts msigdb geneset to fgsea compatible
#' @export
#'
fgsea_msigdb <- function(msigdbgeneset){
  x1 <- msigdbgeneset %>% dplyr::select(gs_name, entrez_gene)
  x1 <- x1 %>% group_by(gs_name) %>% nest()
  geneset <- lapply(x1$data , function(x){as.character(x[[1]])})
  names(geneset) <- x1$gs_name
  return(geneset)
}

#' retrieve several msigdbr collections as fgsea compatible lists
#' @export
fgsea_msigdb_collections <- function(
  msigCollection,
  species = "Homo sapiens"){
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


#' run for all contrasts
#' @export
#' @param allrnk list of rank arrays
#' @param geneSet single list geneset (e.g. )
#' @param nperm number permutations
#' @param minSize minimum geneset size
#' @param maxSize maximum geneset size
run_fgsea_for_allContrasts <- function(allrnk,
                                       geneSet,
                                       nperm = 10000,
                                       minSize = 25,
                                       maxSize = 500){
  fgseaRes <- vector(mode = "list", length = length(allrnk))
  for (i in 1:length(allrnk)) {
    message( paste( names(allrnk)[i],"\n" ) )
    fgseaRes[[i]] <- fgsea::fgsea(pathways = geneSet,
                                  stats    = allrnk[[i]],
                                  nperm = nperm,
                                  minSize  = minSize,
                                  maxSize  = maxSize)

  }

  for (i in 1:length(allrnk)) {
    fgseaRes[[i]]$comparison <- names(allrnk)[i]
  }
  return(fgseaRes)
}

#' used to analyse a single contrast with various genesets
#' @export
#' @param allrnk ranklist
#' @param geneSets list of gene sets
#' @param nperm number of permutations
#' @param minSize minimum geneset size
#' @param maxSize maximum geneset size
#'
run_fgsea_for_allGeneSets <- function(allrnk,
                                      geneSets,
                                      nperm = 10000,
                                      minSize = 25,
                                      maxSize = 500){
  fgseaRes <- vector(mode = "list", length = length(geneSets))
  for (i in 1:length(geneSets)) {
    message( paste( names(geneSets)[i],"\n" ) )
    fgseaResult <- fgsea::fgsea(pathways = geneSets[[i]],
                                stats    = allrnk,
                                nperm = nperm,
                                minSize  = minSize,
                                maxSize  = maxSize)
    if (nrow(fgseaResult) == 0) {
      next()
    }
    fgseaResult$GS <-  names(geneSets)[i]
    relevantResult <- fgseaResult %>%
      dplyr::relocate(nMoreExtreme ,
                      pval,
                      ES,
                      leadingEdge,
                      .after = size)
    fgseaResult <- dplyr::relocate(fgseaResult, GS, .before = pathway)
    fgseaResult <- dplyr::rename(fgseaResult, FDR = padj)
    fgseaRes[[i]] <- fgseaResult

  }
  names(fgseaRes) <- names(geneSets)
  return(fgseaRes)
}

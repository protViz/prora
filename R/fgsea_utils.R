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
    df <- select(df, all_of(c(contrast, ids, score)))
    ldf <- group_by(df, !!sym(contrast)) %>% nest()

    allrnk <- vector(mode = "list", length = length(ldf))
    for (i in 1:length(ldf)) {
      ranks <- fgsea_rank(ldf[[i]], ids = ids, score = score)
      allrnk[[i]] <- ranks
    }
    names(allrnk) <- names(ldf)
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

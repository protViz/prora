#' UpSetR wrapper for sigora results
#'
#' @param sigora_res Object returned by the \code{\link{sigoraWrappR}} function
#' @param GPStable Object returned by \code{\link{makeGPS_wrappR}} function, setting \code{dev=TRUE}
#' @param ... other parameters to \code{upset}
#'
#' @importFrom UpSetR fromList upset
#' @importFrom dplyr inner_join filter select mutate
#' @importFrom tidyr spread
#' @importFrom S4Vectors na.omit
#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @export
sigora_upsetR <- function(sigora_res, GPStable, ...) {
  mergeR <- function(sigora_res, GPStable) {
    tab1 <- sigora_res$sigora$summary_results
    tab2 <- sigora_res$data
    colnames(tab1)[1] <- "pathwayId"
    colnames(tab2) <- c("gene", "fc")
    tab3 <- dplyr::inner_join(GPStable, tab1)
    tab4 <- dplyr::inner_join(tab3, tab2) %>% filter(Bonferroni < 0.05)
    return(tab4)
  }

  df <- mergeR(sigora_res, GPStable)

  if(any(dim(df)==0)) return(NULL)

  toplot <- df %>%
    dplyr::select(pathwayId, gene) %>%
    dplyr::mutate(ID = seq_len(nrow(.data))) %>%
    tidyr::spread(pathwayId, gene) %>%
    dplyr::select(-ID) %>%
    as.list() %>%
    lapply(na.omit)

  if(length(toplot)==1) stop("UpSetR plot cannot be displayed. Only one pathway enriched.")

  UpSetR::fromList(toplot) %>%
    UpSetR::upset(mb.ratio = c(0.7, 0.3), ...=...)
}


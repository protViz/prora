#' UpSetR wrapper for sigora results
#' @export sigora_upsetR sigora_upsetR
#' @import UpSetR

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

  test <- mergeR(sigora_example, GPStab)

  test %>%
    dplyr::select(pathwayId, gene) %>%
    dplyr::mutate(ID = 1:nrow(.)) %>%
    tidyr::spread(pathwayId, gene) %>%
    dplyr::select(-ID) %>%
    as.list() %>%
    lapply(na.omit) %>%
    UpSetR::fromList() %>%
    UpSetR::upset(mb.ratio = c(0.7, 0.3), ...=...)
}


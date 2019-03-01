#' Visualisation of expression patterns in significantly enriched pathways
#'
#' @param sigora_res Object returned by the \code{\link{sigoraWrappR}} function
#' @param GPS_table Object returned by \code{\link{makeGPS_wrappR}} function, setting \code{dev=TRUE}
#'
#' @return Returns a \code{ggplot} object heatmap
#'
#' @examples
#' data("sigora_example", package = "fgczgseaora")
#' data("GPStab", package = "fgczgseaora")
#' sigora_heatmap(sigora_example, GPStab)
#'
#' @export sigora_heatmap sigora_heatmap
#'
#' @import ggplot2

sigora_heatmap <- function(sigora_res, GPStable) {
  tab1 <- sigora_res$sigora$summary_results
  tab2 <- sigora_res$data
  colnames(tab1)[1] <- "pathwayId"
  colnames(tab2) <- c("gene", "fc")
  tab3 <- dplyr::inner_join(GPStable, tab1)
  tab4 <- dplyr::inner_join(tab3, tab2) %>% filter(Bonferroni < 0.05)
  stopifnot(all(dim(tab4)>0))
    ggplot(tab4, aes(x = pathwayId, y = gene, fill = fc)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red") +
    labs(
      title = "Expression Patterns",
      subtitle = "Fold changes of genes from significantly enriched pathways",
      fill = "Fold Change"
    ) +
    xlab("Pathway") +
    ylab("Gene") +
    theme_bw() %>% return()
}


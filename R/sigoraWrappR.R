#' sigora Wrapper function
#'
#'
#' Provides a wrapper for \code{sigora::sigora()} combined with \code{sigora::ora()} and
#' returns a \code{list} used for the generation of an \code{.Rmd} report.
#'
#' @param data input data.frame (at least two columns, first column containing IDs,
#' other columns numerical ranks, i.e. fold changes)
#' @param threshold fold change threshold above which (in absolute terms) a protein
#' is considered differentially regulated
#' @param score_col Name of the fold change column, in case the input file contains
#' multiple contrasts
#' @param GPSrepos GPS repository used as background, can be
#' generated via \code{\link{sigoraWrappR}}
#' @param db database used for the generation of the GPS repository
#' @param greater_than Logical. Whether fc_threshold should be applied as
#' greater_than (default is \code{TRUE}) or strictly less than (\code{FALSE})
#' @param idmap id mapping table, dafault sigora::idmap
#'
#' @return Returns a \code{list} containing the following elements:
#' \itemize{
#'    \item \code{sigora}: \code{data.frame} containing sigora results
#'    \item \code{ora}: \code{data.frame} containing ora results
#'    \item \code{fc_threshold}: fc_threshold used
#'    \item \code{GPS_repository}: GPS repository used, to be reported
#'    \item \code{database}: Database used for generating the GPS repository
#'    \item \code{data}: \code{data.frame} containing Uniprot IDs and fold changes
#'    \item \code{proteinsAfterFiltering}: numeric, number of proteins after fc filtering
#' }
#'
#' @export sigoraWrappR sigoraWrappR
#'
#' @examples
#'
#' library(prora)
#' df <- get_UniprotID_from_fasta_header(prora::exampleContrastData)
#' myGPSrepo <- makeGPS_wrappR(ids = df$UniprotID,target = "KEGG")
#' res <- sigoraWrappR(df,score_col = "estimate", GPSrepos = myGPSrepo$gps,
#'                      threshold = 0.5)
#'
#' myGPSrepoGO <- makeGPS_wrappR(df$UniprotID, target = "GO")
#' res <- sigoraWrappR(df,score_col = "estimate", GPSrepos = myGPSrepoGO$gps,
#'                      threshold = 0.5)
#' \dontrun{
#'   res$ora
#' res$sigora$summary_results
#' res$data
#' GPSrepo <-  myGPSrepoGO$gps
#' }
sigoraWrappR <-
  function(data,
           threshold = 0.5,
           score_col = "",
           GPSrepos = sigora::load_data('kegH'),
           db = "",
           greater_than = TRUE,
           idmap = sigora::idmap) {
    enriched <- NULL
    if (greater_than) {
      enriched <- data[data[, score_col] >= threshold, ]
    } else {
      enriched <- data[data[, score_col] < threshold, ]
    }
    sigora_res <-
      sigora::sigora(
        GPSrepo = GPSrepos,
        level = 2,
        queryList = enriched$UniprotID,
        idmap = idmap
      )
    ora_res <-
      sigora::ora(geneList = enriched$UniprotID,
                  GPSrepo = GPSrepos,
                  idmap = idmap)

    output <- list(
      sigora = sigora_res,
      ora = ora_res,
      threshold = threshold,
      GPSrepository = GPSrepos,
      database = db,
      data = data[, c("UniprotID", score_col)],
      proteinsAfterFiltering = nrow(enriched),
      direction = ifelse(greater_than, yes = "greater than", no = "less than")
    )
    return(output)
  }


.mergeR <- function(sigora_res, GPStable) {
  tab1 <- sigora_res$sigora$summary_results
  tab2 <- sigora_res$data
  colnames(tab1)[1] <- "pathwayId"
  colnames(tab2) <- c("gene", "fc")
  tab3 <- dplyr::inner_join(GPStable, tab1)
  tab4 <-
    dplyr::inner_join(tab3, tab2) %>% filter(!!sym("Bonferroni") < 0.05)
  return(tab4)
}

#' UpSetR wrapper for sigora results
#'
#' @param sigora_res Object returned by the \code{\link{sigoraWrappR}} function
#' @param GPStable Object returned by \code{\link{makeGPS_wrappR}} function, setting \code{dev=TRUE}
#' @param ... other parameters to \code{upset}
#'
#' @export
#' @examples
#' sigora_res <- prora::sigora_example
#' GPStable <- prora::GPStab
#'
#' sigora_upsetR(sigora_res, GPStable)
#'
#'
#' df <- get_UniprotID_from_fasta_header(prora::exampleContrastData)
#' myGPSrepo <- makeGPS_wrappR(ids = df$UniprotID)
#' names(myGPSrepo)
#' res <- sigoraWrappR(df,score_col = "estimate", GPSrepos = myGPSrepo$gps,
#'                      threshold = 0.5)
#' sigora_upsetR(res, myGPSrepo$gpsTable)
#'
sigora_upsetR <- function(sigora_res, GPStable, ...) {
  df <- .mergeR(sigora_res, GPStable)
  if (any(dim(df) == 0))
    return(NULL)
  toplot <- df %>%
    dplyr::select(!!sym("pathwayId"),!!sym("gene")) %>%
    dplyr::mutate(ID = seq_len(nrow(df))) %>%
    tidyr::spread(!!sym("pathwayId"),!!sym("gene")) %>%
    dplyr::select(-!!sym("ID")) %>%
    as.list() %>%
    lapply(na.omit)

  if (length(toplot) == 1)
    stop("UpSetR plot cannot be displayed. Only one pathway enriched.")
  UpSetR::upset(UpSetR::fromList(toplot), mb.ratio = c(0.7, 0.3), ... = ...)
}

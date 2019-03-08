#' sigora Wrapper function
#'
#'
#' Provides a wrapper for \code{sigora::sigora()} combined with \code{sigora::ora()} and returns a \code{list} used for the generation of an \code{.Rmd} report.
#'
#' @param fc_threshold fold change threshold above which (in absolute terms) a protein is considered differentially regulated
#' @param fc_col Name of the fold change column, in case the input file contains multiple contrasts
#' @param GPSrepos GPS repository used as background, can be generated via \code{\link{sigoraWrappR}}
#' @param df data.frame input (at least two columns, first column being gene symbols, other columns being fold changes)
#' @param db database used for the generation of the GPS repository
#'
#' @return Returns a \code{list} containing the following elements:
#' \itemize{
#'    \item \code{sigora}: \code{data.frame} containing sigora results
#'    \item \code{ora}: \code{data.frame} containing ora results
#'    \item \code{fc_threshold}: fc_threshold used
#'    \item \code{GPS_repository}: GPS repository used, to be reported
#'    \item \code{database}: Database used for generating the GPS repository
#'    \item \code{data}: \code{data.frame} containing gene symbols and fold changes
#'    \item \code{proteinsAfterFiltering}: numeric, number of proteins after fc filtering
#' }
#'
#' @examples
#' data("exampleContrastData", package = "fgczgseaora")
#' data("idmap", package = "sigora")
#' df <- getSymbolFromSwissprotID(exampleContrastData)
#' sigoraWrappR(fc_col = colnames(df)[7], GPSrepos = sigora::kegH, df = df)
#'
#' @import tidyverse sigora org.Hs.eg.db GO.db reactome.db dplyr AnnotationDbi tidyr
#' @export sigoraWrappR sigoraWrappR

sigoraWrappR <-
  function(fc_threshold = 0.5,
           fc_col = "",
           GPSrepos = sigora::kegH,
           df,
           db = "") {
    if(fc_threshold >= 0) {
      enriched <- df[df[, fc_col] >= fc_threshold,]
    } else {
      enriched <- df[df[, fc_col] <= fc_threshold,]
    }
    sigora_res <-
      sigora::sigora(GPSrepo = GPSrepos,
             level = 5,
             queryList = enriched$Symbol)
    ora_res <-
      sigora::ora(geneList = enriched$Symbol,
          GPSrepo = GPSrepos)
    output <- list(
      sigora = sigora_res,
      ora = ora_res,
      fc_threshold = fc_threshold,
      GPSrepository = GPSrepos,
      database = db,
      data = df[, c("Symbol", fc_col)],
      proteinsAfterFiltering = nrow(df[df[, fc_col] >= fc_threshold,])
    )
    return(output)
  }

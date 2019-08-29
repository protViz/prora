#' sigora Wrapper function
#'
#'
#' Provides a wrapper for \code{sigora::sigora()} combined with \code{sigora::ora()} and
#' returns a \code{list} used for the generation of an \code{.Rmd} report.
#'
#' @param fc_threshold fold change threshold above which (in absolute terms) a protein
#' is considered differentially regulated
#' @param fc_col Name of the fold change column, in case the input file contains
#' multiple contrasts
#' @param GPSrepos GPS repository used as background, can be
#' generated via \code{\link{sigoraWrappR}}
#' @param df input data.frame (at least two columns, first column containing IDs,
#' other columns numerical ranks, i.e. fold changes)
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
#' @examples
#' data("exampleContrastData", package = "fgczgseaora")
#' data("idmap", package = "sigora")
#' df <- getUniprotFromFastaHeader(exampleContrastData)
#' myGPSrepo <- makeGPS_wrappR(ids = df$UniprotID)
#' res <- sigoraWrappR(fc_col = "estimate", GPSrepos = myGPSrepo,
#'                     df = df, fc_threshold = 0.5)
#'
#' @importFrom sigora sigora ora
#' @export sigoraWrappR sigoraWrappR

sigoraWrappR <-
  function(fc_threshold = 0.5,
           fc_col = "",
           GPSrepos = sigora::kegH,
           df,
           db = "",
           greater_than = TRUE,
           idmap = sigora::idmap) {
    if(greater_than) {
      enriched <- df[df[, fc_col] >= fc_threshold,]
    } else {
      enriched <- df[df[, fc_col] < fc_threshold,]
    }
    sigora_res <-
      sigora::sigora(GPSrepo = GPSrepos,
             level = 2,
             queryList = enriched$UniprotID,
             idmap = idmap)
    ora_res <-
      sigora::ora(geneList = enriched$UniprotID,
          GPSrepo = GPSrepos)
    output <- list(
      sigora = sigora_res,
      ora = ora_res,
      fc_threshold = fc_threshold,
      GPSrepository = GPSrepos,
      database = db,
      data = df[, c("UniprotID", fc_col)],
      proteinsAfterFiltering = nrow(df[df[, fc_col] >= fc_threshold, ]),
      direction = ifelse(greater_than, yes = "greater than", no = "less than")
    )
    return(output)
  }

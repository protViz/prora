#' Contrast effect estimates
#'
#' Example data set containing effect estimates and p values from different contrasts.
#'
#' @docType data
#' @name exampleContrastData
#' @usage data(exampleContrastData)
NULL

#' Examplary Uniprot IDs
#'
#' Gene symbols for running examples to produce GPS repositories using \code{\link{makeGPS_wrappR}}
#'
#' @docType data
#' @name exampleUniprotIDs
#' @usage data(exampleUniprotIDs)
NULL

#' Exemplary GPS table
#'
#' GPS table linking GOID, pathway name and genes
#' \itemize{
#'   \item \code{pathwayId.} Pathway identifier, in this case a gene ontology id
#'   \item \code{pathwayName.} Name of the pathway
#'   \item \code{gene.} Gene symbol
#' }
#'
#' @docType data
#' @name GPStab
#' @usage data(GPStab)
NULL

#' Examplary GPS repository
#' @docType data
#' @name myGPSrepo
#' @usage data(myGPSrepo)
NULL

#' Examplary GPS repository
#' @docType data
#' @name myGPSrepo
#' @usage data(myGPSrepo)
NULL


#' Exemplary GSEA output
#'
#' GSEA output from the \code{run_GSEA.R} script provided in the \code{inst/RunScripts/} directory
#'
#' \itemize{
#'   \item \code{organism.} The organism used in the analysis
#'   \item \code{target.} Target database
#'   \item \code{outputDir.} Output directory of where the results are stores
#'   \item \code{input_data.} Input data containing gene symbols and effect estimates
#'   \item \code{GSEA_res.} Results from \code{{WebGestaltR::WebGestaltR}}
#' }
#'
#' @docType data
#' @name GSEA
#' @usage data(GSEA)
NULL

#' examplary GSEA output
#' @docType data
#' @name GSEAex1
#' @usage data(GSEAex1)
NULL

#' Exemplary sigora output
#'
#' sigora output from the \code{\link{sigoraWrappR}} function containing the following entries
#'
#' \itemize{
#'   \item \code{sigora.} Result table from \code{{sigora::sigora}}
#'   \item \code{ora.} Result table from \code{{sigora::ora}}
#'   \item \code{fc_threshold.} fold change threshold used for the analysis
#'   \item \code{GPSrepository.} GPS repository used for the analysis, generated by \code{\link{makeGPS_wrappR}}
#'   \item \code{database.} Target database
#'   \item \code{data.} Original input data
#'   \item \code{proteinsAfterFiltering.} The number of proteins meeting the fold change threshold
#' }
#'
#' @docType data
#' @name sigora_example
#' @usage data(sigora_example)
NULL


#' workflow for sigORA
#' @param data data
#' @param ID_col column name containing IDs
#' @param score_col column name containing estimates
#' @param fc_threshold threshold for estimate
#' @param greater flag whether to filter > threshold or < threshold
#' @param target target database, default: "GO"
#' @param outdir output directory
#' @export
#'
#' @examples
#'
#' library(prora)
#' fc_estimates <- prora::exampleContrastData
#' filtered_dd <- get_UniprotID_from_fasta_header(fc_estimates, idcolumn = "protein_Id")
#' myGPSrepo <- makeGPS_wrappR(filtered_dd$UniprotID, target = "GO")
#'
#' #usethis::use_data(myGPSrepo,overwrite = TRUE)
#' undebug(runSIGORA)
#' undebug(sigoraWrappR)
#' res <- runSIGORA(filtered_dd, myGPSrepo = prora::myGPSrepo, DEBUG=TRUE)
#' #rmarkdown::render(res$rmarkdownPath,bookdown::html_document2(number_sections = FALSE),params = res$sigoraData,clean = TRUE)
#'
runSIGORA <-
  function(data,
           target = "GO",
           score_col = "estimate",
           ID_col = "UniprotID",
           threshold = 0.5,
           outdir = "sigORA",
           greater = TRUE,
           myGPSrepo = NULL,
           DEBUG=FALSE
           ) {

    if(is.null(myGPSrepo)){
      myGPSrepo <-
        makeGPS_wrappR(data[[ID_col]], target = target)
    }

    sigora_res <-
      try(sigoraWrappR(
        data,
        threshold = threshold,
        score_col = score_col,
        GPSrepos = myGPSrepo$gps,
        db = target,
        greater_than = greater
      ))

    if (inherits(sigora_res, "try-error")) {
      message("No enriched pathways. No report will be generated.")
      return(invisible(NULL))
    }

    sigoraData <- list(
      results = sigora_res,
      GPStable = myGPSrepo$gpsTable ,
      direction_greater = greater
    )

    outdir <- paste0(outdir, "_", target)
    rmarkdownPath <- file.path(outdir, "sigora.Rmd")

    if (!DEBUG){
      if (!dir.exists(outdir)){
        dir.create(outdir)
      }

      rmarkdownPath_src <- file.path(
        find.package("prora"),
        "rmarkdown_reports/sigora.Rmd"
      )
      if(!file.copy(rmarkdownPath_src,
                    rmarkdownPath,
                    overwrite = TRUE
      )){
        warning("could not copy",rmarkdownPath_src, "to ", rmarkdownPath)
      }

      bibpath <- file.path(outdir, "bibliography.bib")
      bibpath_src <- file.path(
        find.package("prora"),
        "rmarkdown_reports/bibliography.bib"
      )
      if(!file.copy(
        bibpath_src,
        bibpath,
        overwrite = TRUE
      )){
        warning("could not copy",bibpath_src, "to ", bibpath)
      }
      rmarkdown::render(
        rmarkdownPath,
        bookdown::html_document2(number_sections = FALSE),
        params = sigoraData,
        clean = TRUE
      )
    }


    return(list(sigoraData = sigoraData, rmarkdownPath =  rmarkdownPath))
  }

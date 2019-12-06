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
#' library(fgczgseaora)
#' fc_estimates <- fgczgseaora::exampleContrastData
#' filtered_dd <- get_UniprotID_from_fasta_header(fc_estimates, idcolumn = "protein_Id")
#'
#' myGPSrepo <- makeGPS_wrappR(filtered_dd$UniprotID, target = "GO")
#' usethis::use_data(myGPSrepo)
#' names(myGPSrepo)
#' res <- runSIGORA(filtered_dd, myGPSrepo = myGPSrepo)
#' names(res)
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
           myGPSrepo = NULL
           ) {
    outdir <- paste0(outdir, "_", target)

    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }

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

    if(!dir.exists(outdir)){
      dir.create(outdir)
    }
    rmarkdownPath <- file.path(outdir, "sigora.Rmd")

    file.copy(
      file.path(
        find.package("fgczgseaora"),
        "rmarkdown_reports/sigora.Rmd"
      ),
      rmarkdownPath,
      overwrite = TRUE
    )

    bibpath <- file.path(outdir, "bibliography.bib")
    file.copy(
      file.path(
        find.package("fgczgseaora"),
        "rmarkdown_reports/bibliography.bib"
      ),
      bibpath,
      overwrite = TRUE
    )

    sigoraData <- list(
      results = sigora_res,
      GPStable = myGPSrepo$gpsTable ,
      direction_greater = greater
    )

    rmarkdown::render(
      rmarkdownPath,
      bookdown::html_document2(number_sections = FALSE),
      params = sigoraData,
      clean = TRUE
    )
    return(list(sigoraData = sigoraData, rmarkdownPath =  rmarkdownPath))
  }

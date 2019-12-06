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
#' fc_estimates <- fgczgseaora::exampleContrastData
#' filtered_dd <- get_UniprotID_from_fasta_header(fc_estimates, idcolumn = "protein_Id")
#'
#' runSIGORA(filtered_dd)
runSIGORA <-
  function(data,
           target = "GO",
           score_col = "estimate",
           ID_col = "UniprotID",
           threshold = 0.5,
           outdir = "sigORA",
           greater = TRUE) {
    outdir <- paste0(outdir, "_", target)

    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }

    myGPSrepo <-
      makeGPS_wrappR(data[[ID_col]], target = target)

    sigora_res <-
      try(sigoraWrappR(
        data,
        threshold = threshold,
        score_col = score_col,
        GPSrepos = myGPSrepo,
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
      GPStable = GPStab,
      direction_greater = greater
    )
    rmarkdown::render(
      rmarkdownPath,
      bookdown::html_document2(number_sections = FALSE),
      params = ,
      clean = TRUE
    )
    return(sigoraData)
  }

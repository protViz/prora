# Helper functions for run scripts

#' @importFrom rlang sym
#' @importFrom dplyr filter
.apply_threshold <-
  function(df, th, col = "Score", greater = TRUE) {
    if (greater) {
      df %>%
        filter(!!sym(col) >= th) -> out
    } else {
      df %>%
        filter(!!sym(col) <= th) -> out
    }
    return(out)
  }

#' workflow for GSEA
#'
#' @param data data
#' @param fpath file path to write to
#' @param ID_col column name containing IDs
#' @param fc_col column name containing estimates
#' @param organism organism
#' @param target target database, default: geneontology_Biological_Process
#' @param nperm number of permutations
#' @param outdir output directory
#' @param interestGeneType what type of identifier default : "uniprotswissprot"
#' @param contrast_name for pretty printing
#' @importFrom WebGestaltR WebGestaltR
#' @importFrom readr read_delim
#' @importFrom tidyr separate_rows
#' @export
#'
runGSEA <- function(data,
                    fpath,
                    ID_col = "UniprotID",
                    fc_col = "estimate",
                    organism = "hsapiens",
                    target = "geneontology_Biological_Process",
                    nperm = 10,
                    outdir = "GSEA",
                    interestGeneType = "uniprotswissprot",
                    contrast_name = fpath)
{
  outdir <- file.path(outdir, target)

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  ranktable <- data %>%
    dplyr::select(!!sym(ID_col), Score = !!sym(fc_col))

  GSEA_res <-
    WebGestaltR(
      enrichMethod = "GSEA",
      organism = organism,
      enrichDatabase = target,
      interestGene = ranktable,
      interestGeneType = interestGeneType,
      outputDirectory = outdir,
      isOutput = TRUE,
      perNum = nperm,
      projectName = fpath
    )

  if(is.null(GSEA_res)){
    warning("!!! no results returned by WebGestaltR !!!", GSEA_res)
    return(NULL)
  }

  rdataPath <- file.path(outdir, paste0("Project_", fpath), "GSEA_res.Rdata")

  message("storing GSEA_res to: ", rdataPath)
  saveRDS(GSEA_res, file = rdataPath)

  f_mappingTable <- file.path(
    outdir,
    paste0("Project_", fpath),
    paste0("interestingID_mappingTable_", fpath, ".txt")
  )

  mappingTable <- read_delim(f_mappingTable, delim = "\t") %>%
    mutate(entrezgene = as.character(!!sym("entrezgene")))
  GSEA_res_sep <-  GSEA_res %>%
    separate_rows(!!sym("leadingEdgeId"), sep = ";")
  merged_data <- inner_join(mappingTable,
                            GSEA_res_sep,
                            by = c("entrezgene" = "leadingEdgeId"))

  readr::write_delim(
    merged_data,
    path = file.path(outdir, paste0("Project_", fpath), "merged_data.tsv"),
    delim = "\t"
  )

  GSEA <- list(
    outdir = outdir,
    target = target,
    fpath = fpath,
    organism = organism,
    input_data = ranktable,
    merged_data = merged_data,
    nperm = nperm,
    interestGeneType = interestGeneType,
    contrast_name = contrast_name
  )

  rdataPath <- file.path(outdir, paste0("Project_", fpath), "GSEA.Rdata")
  message("storing GSEA.rdata to: ", rdataPath)
  saveRDS(GSEA, file= rdataPath)

  rmarkdownPath <-
    file.path(outdir, paste0("Project_", fpath), "GSEA.Rmd")
  bibpath <-
    file.path(outdir, paste0("Project_", fpath), "bibliography.bib")


  file.copy(
    file.path(find.package("fgczgseaora"), "rmarkdown_reports/GSEA.Rmd"),
    rmarkdownPath,
    overwrite = TRUE
  )

  file.copy(
    file.path(
      find.package("fgczgseaora"),
      "rmarkdown_reports/bibliography.bib"
    ),
    bibpath,
    overwrite = TRUE
  )

  message("rendering :", rmarkdownPath)
  rmarkdown::render(
    rmarkdownPath,
    bookdown::html_document2(number_sections = FALSE),
    params = list(GSEA = GSEA),
    clean = TRUE,
    envir = new.env()
  )
  return(GSEA)
}

#' workflow for sigORA
#' @param data data
#' @param ID_col column name containing IDs
#' @param fc_col column name containing estimates
#' @param fc_threshold threshold for estimate
#' @param greater flag whether to filter > threshold or < threshold
#' @param target target database, default: "GO"
#' @param outdir output directory
#' @export
runSIGORA <-
  function(data,
           target = "GO",
           fc_col = "estimate",
           ID_col = "UniprotID",
           fc_threshold = 0.5,
           outdir = "sigORA",
           greater = TRUE) {
    outdir <- paste0(outdir, "_", target)

    if (!dir.exists(outdir)) {
      dir.create(outdir, recursive = TRUE)
    }

    GPStab <-
      makeGPS_wrappR(data[[ID_col]], target = target, dev = TRUE)

    myGPSrepo <-
      makeGPS_wrappR(data[[ID_col]], target = target)

    sigora_res <-
      try(sigoraWrappR(
        fc_threshold = fc_threshold,
        fc_col = fc_col,
        df = data,
        GPSrepos = myGPSrepo,
        db = target,
        greater_than = greater
      ))

    if (inherits(sigora_res, "try-error")) {
      message("No enriched pathways. No report will be generated.")
      return(invisible(NULL))
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

#' workflow for ORA
#' @param data data
#' @param fpath file path to write to
#' @param ID_col column name containing IDs
#' @param fc_col column name containing estimates
#' @param organism organism
#' @param target target database, default: "geneontology_Biological_Process"
#' @param threshold threshold for estimate
#' @param nperm number of permutations
#' @param outdir output directory
#' @param greater indicating direction of threshold
#' @param interestGeneType what type of identifier default : "uniprotswissprot"
#'
#' @importFrom WebGestaltR WebGestaltR
#' @export
runWebGestaltORA <- function(data,
                             fpath,
                             organism = "hsapiens",
                             ID_col = "UniprotID",
                             target = "geneontology_Biological_Process",
                             threshold = 0.5,
                             greater = TRUE,
                             nperm = 10,
                             fc_col = "estimate",
                             outdir = "WebGestalt_ORA",
                             interestGeneType = "uniprotswissprot",
                             contrast_name = fpath) {
  outdir <- file.path(outdir, target)

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  dat <- data %>%
    dplyr::select(!!sym(ID_col), Score = !!sym(fc_col))

  ranktable <-
    .apply_threshold(df = dat,
                     th = threshold,
                     greater = greater)

  message("There are ", nrow(dat), " proteins in the background\n\n")
  message("There are ",nrow(ranktable) , " proteins in the selected set\n\n" )

  if(nrow(ranktable) == 0){
    message("NO proteins in the subset!\n")
    return(0)
  }

  ORA_res <- tryCatch(
    WebGestaltR(
      enrichMethod = "ORA",
      organism = organism,
      enrichDatabase = target,
      interestGene = ranktable[[ID_col]],
      referenceGene = data[[ID_col]],
      interestGeneType = interestGeneType,
      referenceGeneType = interestGeneType,
      outputDirectory = outdir,
      isOutput = TRUE,
      perNum = nperm,
      projectName = fpath
    ), error = function(e){message(e) ; return(NULL) })
  message("\n\n Finished ORA \n\n")

  ORA <- list(
    outdir = outdir,
    target = target,
    fpath = fpath,
    organism = organism,
    input_data = data,
    threshold = threshold,
    greater = greater,
    merged_data = ORA_res,
    nperm = nperm,
    interestGeneType = interestGeneType,
    contrast_name = contrast_name
  )

  return(ORA_res)
}

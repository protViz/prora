# Helper functions for run scripts

.apply_threshold <-
  function(df,
           th,
           greater = TRUE) {
    if (greater) {
      df %>%
        dplyr::filter(Score >= th) -> out
    } else {
      df %>%
        dplyr::filter(Score <= th) -> out
    }
    return(out)
  }

.runGSEA <- function(data,
                    fpath,
                    ID_col = "UniprotID",
                    fc_col = "estimate",
                    organism = "hsapiens",
                    target = "geneontology_Biological_Process",
                    nperm = 10,
                    outdir = "GSEA") {
  outdir <- file.path(outdir, target)

  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }

  ranktable <- data %>%
    dplyr::select(!!sym(ID_col), Score = !!sym(fc_col))

  GSEA_res <-
    WebGestaltR(
      enrichMethod = "GSEA",
      organism = organism,
      enrichDatabase = target,
      interestGene = ranktable,
      interestGeneType = "uniprotswissprot",
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
  saveRDS(GSEA_res, file= rdataPath)

  f_mappingTable <- file.path(
    outdir,
    paste0("Project_", fpath),
    paste0("interestingID_mappingTable_", fpath, ".txt")
  )
  mappingTable <-
    read_delim(f_mappingTable,
               delim = "\t")

  mappingTable %>% mutate(entrezgene = as.character(entrezgene)) -> mappingTable
  GSEA_res_sep <-
    GSEA_res %>% separate_rows(leadingEdgeId, sep = ";")
  merged_data <- inner_join(mappingTable,
                            GSEA_res_sep,
                            by = c("entrezgene" = "leadingEdgeId"))

  readr::write_delim(
    merged_data,
    path = file.path(outdir, paste0("Project_", fpath), "merged_data.tsv"),
    delim = "\t"
  )

  GSEA <- list(
    organism = organism,
    target = target,
    input_data = ranktable,
    output_dir = fpath,
    merged_data = merged_data,
    nperm = nperm
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
    clean = TRUE
  )
  return(GSEA)
}

.runSIGORA <-
  function(data,
           target = "GO",
           fc_col = "estimate",
           ID_col = "UniprotID",
           fc_threshold = 0.5,
           outdir = "sigORA",
           greater = TRUE) {
    outdir <- paste0(outdir, "_", target)

    if (!dir.exists(outdir)) {
      dir.create(outdir)
    }

    GPStab <-
      makeGPS_wrappR(data[[ID_col]], target = target, dev = TRUE)

    myGPSrepo <-
      makeGPS_wrappR(data[[ID_col]], target = target)

    sigora_res <-
      sigoraWrappR(
        fc_threshold = fc_threshold,
        fc_col = fc_col,
        df = data,
        GPSrepos = myGPSrepo,
        db = target,
        greater_than = greater
      )


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


.runWebGestaltORA <- function(data,
                             fpath,
                             organism = "hsapiens",
                             ID_col = "UniprotID",
                             target = "geneontology_Biological_Process",
                             threshold = 0.5,
                             greater = TRUE,
                             nperm = 10,
                             fc_col = "estimate",
                             outdir = "WebGestalt_ORA") {
  outdir <- file.path(outdir, target)

  if (!dir.exists(outdir)) {
    dir.create(outdir)
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

  ORA_res <- try(
    WebGestaltR(
      enrichMethod = "ORA",
      organism = organism,
      enrichDatabase = target,
      interestGene = ranktable[[ID_col]],
      referenceGene = data[[ID_col]],
      interestGeneType = "uniprotswissprot",
      referenceGeneType = "uniprotswissprot",
      outputDirectory = outdir,
      isOutput = TRUE,
      perNum = nperm,
      projectName = fpath
    ))
  return(ORA_res)
}

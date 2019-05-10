#' @export runGSEAlong
runGSEAlong <- function(contrast,
                        organism = "hsapiens",
                        ID_col = "UniprotID",
                        target = "geneontology_Biological_Process",
                        map_col = "GO",
                        nperm = 10,
                        fc_col = "estimate",
                        contrast_col = "lhs",
                        outdir = "GSEA") {
  outdir <- paste0(outdir, "_", target)
  fpath <- make.names(contrast)

  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }

  ranktable <- ddd %>%
    dplyr::filter(!!sym(contrast_col) == contrast) %>%
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

  rmarkdown::render(
    rmarkdownPath,
    bookdown::html_document2(number_sections = FALSE),
    params = list(GSEA = GSEA),
    clean = TRUE
  )
}

#' @export runSIGORAlong
runSIGORAlong <-
  function(contrast,
           target = "GO",
           fc_col = "estimate",
           fc_threshold = 0.5,
           outdir = "sigORA",
           contrast_col = con_col,
           greater = TRUE) {
    outdir <- paste0(outdir, "_", target)
    fpath <- make.names(contrast)

    if (!dir.exists(outdir)) {
      dir.create(outdir)
    }

    if (!dir.exists(file.path(outdir, fpath))) {
      dir.create(file.path(outdir, fpath))
    }

    dat <- ddd %>%
      dplyr::filter(!!sym(contrast_col) == contrast)

    GPStab <-
      makeGPS_wrappR(dat$UniprotID, target = target, dev = TRUE)

    myGPSrepo <-
      makeGPS_wrappR(dat$UniprotID, target = target)

    sigora_res <-
      sigoraWrappR(
        fc_threshold = fc_threshold,
        fc_col = fc_col,
        df = dat,
        GPSrepos = myGPSrepo,
        db = target,
        greater_than = greater
      )

    p1 <- try(sigora_heatmap(sigora_example, GPStab))

    rmarkdownPath <- file.path(outdir, fpath, "sigora.Rmd")

    bibpath <- file.path(outdir, fpath, "bibliography.bib")

    file.copy(
      file.path(
        find.package("fgczgseaora"),
        "rmarkdown_reports/sigora.Rmd"
      ),
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

    rmarkdown::render(
      rmarkdownPath,
      bookdown::html_document2(number_sections = FALSE),
      params = list(
        results = sigora_res,
        plot1 = p1,
        GPStable = GPStab,
        direction_greater = greater
      ),
      clean = TRUE
    )
  }

#' @export apply_threshold
apply_threshold <-
  function(df,
           th,
           greater = TRUE) {
    if (greater) {
      df %>%
        dplyr::filter(Score > th) -> out
    } else {
      df %>%
        dplyr::filter(Score <= th) -> out
    }
    return(out)
  }

#' @export runWebGestaltORAlong
runWebGestaltORAlong <- function(contrast,
                                 organism = "hsapiens",
                                 ID_col = "UniprotID",
                                 target = "geneontology_Biological_Process",
                                 map_col = "GO",
                                 threshold = 0.5,
                                 greater = TRUE,
                                 nperm = 10,
                                 contrast_col = con_col,
                                 fc_col = "estimate",
                                 outdir = "WebGestalt_ORA") {
  outdir <- paste0(outdir, "_", target)
  fpath <- make.names(contrast)

  if (!dir.exists(outdir)) {
    dir.create(outdir)
  }

  dat <- ddd %>%
    dplyr::filter(!!sym(contrast_col) == contrast) %>%
    dplyr::select(!!sym(ID_col), Score = !!sym(fc_col))

  ranktable <-
    apply_threshold(df = dat,
                    th = threshold,
                    greater = greater)

  ORA_res <-
    WebGestaltR(
      enrichMethod = "ORA",
      organism = organism,
      enrichDatabase = target,
      interestGene = ranktable$UniprotID,
      referenceGene = ddd$UniprotID,
      interestGeneType = "uniprotswissprot",
      referenceGeneType = "uniprotswissprot",
      outputDirectory = outdir,
      isOutput = TRUE,
      perNum = nperm,
      projectName = fpath
    )
}

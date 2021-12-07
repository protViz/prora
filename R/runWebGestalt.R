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
#' @param score_col column name containing estimates
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
runWebGestaltGSEA <- function(data,
                    fpath,
                    ID_col = "UniprotID",
                    score_col = "estimate",
                    organism = "hsapiens",
                    target = "geneontology_Biological_Process",
                    nperm = 10,
                    outdir = "GSEA",
                    interestGeneType = "uniprotswissprot",
                    contrast_name = fpath)
{
  outdir <- file.path(outdir, target)

  fpath <- gsub("[[:punct:]]", "_", fpath)

  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE)
  }

  ranktable <- data %>%
    dplyr::select(!!sym(ID_col), Score = !!sym(score_col))


  GSEA_res <- WebGestaltR::WebGestaltR(
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

  rankTable <- file.path(outdir, paste0("Project_", fpath), paste0("rankTable.rnk"))
  message("storing ranktable to file: ")
  readr::write_tsv(ranktable, file = rankTable, col_names = FALSE)

  if (is.null(GSEA_res)) {
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
    file.path(find.package("prora"), "rmarkdown_reports/GSEA.Rmd"),
    rmarkdownPath,
    overwrite = TRUE
  )

  file.copy(
    file.path(
      find.package("prora"),
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


#' workflow for ORA
#' @param data data
#' @param fpath file path to write to
#' @param ID_col column name containing IDs
#' @param score_col column name containing estimates
#' @param organism organism
#' @param target target database, default: "geneontology_Biological_Process"
#' @param threshold threshold for estimate
#' @param nperm number of permutations
#' @param outdir output directory
#' @param greater indicating direction of threshold
#' @param interestGeneType what type of identifier default : "uniprotswissprot"
#' @param subdir_name subdirectory name
#' @param contrast_name default fpath
#'
#' @importFrom WebGestaltR WebGestaltR
#' @export
runWebGestaltORA <- function(data,
                             fpath,
                             organism = "hsapiens",
                             ID_col = "UniprotID",
                             target = "geneontology_Biological_Process",
                             threshold = 0.59,
                             greater = TRUE,
                             nperm = 10,
                             score_col = "estimate",
                             outdir = "WebGestalt_ORA",
                             subdir_name = NULL,
                             interestGeneType = "uniprotswissprot",
                             contrast_name = fpath) {

  if (is.null(subdir_name)) {
    subdir_name <- paste0("fc_threshold_",abs(threshold),"_is_greater_",greater)
  }
  subdir_path <- file.path(outdir, subdir_name)

  if (!dir.exists(subdir_path)) {
    message("created directory : ", subdir_path, "\n\n")
    dir.create(subdir_path)
  }

  outdir_path <- file.path(subdir_path, target)
  if (!dir.exists(outdir_path)) {
    message("created directory : ", outdir_path, "\n\n")
    dir.create(outdir_path, recursive = TRUE)
  }

  dat <- data %>%
    dplyr::select(!!sym(ID_col), Score = !!sym(score_col))

  ranktable <-
    .apply_threshold(df = dat,
                     th = threshold,
                     greater = greater)

  message("There are ", nrow(dat), " proteins in the background\n\n")
  message("There are ", nrow(ranktable) , " proteins in the selected set\n\n" )

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
      outputDirectory = outdir_path,
      subdir = subdir_name,
      isOutput = TRUE,
      perNum = nperm,
      projectName = fpath
    ), error = function(e){message(e) ; return(NULL) })
  message("\n\n Finished ORA \n\n")

  ORA <- list(
    outdir = outdir,
    subdir_name = subdir_name,
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

  return(ORA)
}

rm(list = ls())
library(WebGestaltR)
library(tidyverse)
library(fgczgseaora)
library(org.Hs.eg.db)
library(conflicted)


# Setup -------------------------------------------------------------------

# fpath <- "inst/example_data/Contrasts_SignificanceValues_f_Cells_Treatment.csv"
# dd <- read_csv(fpath)

dd <- contrast_data_example

colnames(dd) <- make.names(colnames(dd))

ddd <- getUniprotFromFastaHeader(dd)

organism <- "hsapiens"
ID_col <- "UniprotID"
target <- "geneontology_Biological_Process"
map_col <- "GO"
nperm <- 10
fc_col <- "estimate"
contrast_col <- "lhs"

contrs <- ddd %>%
  distinct(!!sym(contrast_col)) %>%
  pull()


# Function ----------------------------------------------------------------

runGSEAlong <- function(contrast,
                        organism = "hsapiens",
                        ID_col = "UniprotID",
                        target = "geneontology_Biological_Process",
                        map_col = "GO",
                        nperm = 10,
                        fc_col = "estimate",
                        contrast_col = "lhs",
                        outdir = "GSEA") {
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


# Run ---------------------------------------------------------------------

sapply(contrs, runGSEAlong)

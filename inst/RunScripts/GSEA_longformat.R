rm(list = ls())
library(WebGestaltR)
library(tidyverse)
library(fgczgseaora)
library(org.Hs.eg.db)
library(conflicted)

fpath <- "inst/example_data/Contrasts_SignificanceValues_f_Cells_Treatment.csv"

dd <- read_csv(fpath)
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

for (this.contrast in contrs) {
  this.contrast <- contrs

  fpath <- make.names(this.contrast)

  if(!dir.exists(target)){
    dir.create(target)
  }

  ranktable <- ddd %>%
    dplyr::filter(!!sym(contrast_col) == this.contrast) %>%
    dplyr::select(!!sym(ID_col), Score = !!sym(fc_col))

  GSEA_res <-
    WebGestaltR(
      enrichMethod = "GSEA",
      organism = organism,
      enrichDatabase = target,
      interestGene = ranktable,
      interestGeneType = "uniprotswissprot",
      outputDirectory = target,
      isOutput = TRUE,
      perNum = nperm,
      projectName = fpath
    )

  f_mappingTable <- file.path(target,
                              paste0("Project_",fpath),
                              paste0("interestingID_mappingTable_", fpath,".txt"))
  mappingTable <-
    read_delim(f_mappingTable,
               delim = "\t")

  mappingTable %>% mutate(entrezgene = as.character(entrezgene)) -> mappingTable
  GSEA_res_sep <- GSEA_res %>% separate_rows(leadingEdgeId, sep=";")
  merged_data <- inner_join(mappingTable, GSEA_res_sep,
                            by=c("entrezgene" = "leadingEdgeId"))

  GSEA <- list(
    organism = organism,
    target = target,
    input_data = ranktable,
    output_dir = fpath,
    merged_data = merged_data,
    nperm = nperm
  )

  rmarkdownPath <- file.path(target, paste0("Project_",fpath), "GSEA.Rmd")
  bibpath <- file.path(target,paste0("Project_",fpath), "bibliography.bib")


  file.copy(
    file.path(find.package("fgczgseaora"), "rmarkdown_reports/GSEA.Rmd"),
    rmarkdownPath,
    overwrite = TRUE
  )

  file.copy(
    file.path(find.package("fgczgseaora"), "rmarkdown_reports/bibliography.bib"),
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

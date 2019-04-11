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
threshold <- 0.5
direction <- "greater.than"
nperm <- 10

contrast_col <- "lhs"
fc_col <- "estimate"

contrs <- ddd %>%
  distinct(!!sym(contrast_col)) %>%
  pull()

apply_threshold <-
  function(df,
           th,
           alt = c("greater.than", "less.than")) {
    alt <- match.arg(alt)
    if (alt == "greater.than") {
      df %>%
        dplyr::filter(Score > th) -> out
    } else if (alt == "less.than") {
      df %>%
        dplyr::filter(Score <= th) -> out
    }
    return(out)
  }

for (this.contrast in contrs) {

  if (!dir.exists(this.contrast)) {
    dir.create(this.contrast)
  }

  dat <- ddd %>%
    dplyr::filter(!!sym(contrast_col) == this.contrast) %>%
    dplyr::select(!!sym(ID_col), Score = !!sym(fc_col))

  ranktable <-
    apply_threshold(
      df = dat,
      th = threshold,
      alt = direction
    )

  ORA_res <-
    WebGestaltR(
      enrichMethod = "ORA",
      organism = organism,
      enrichDatabase = target,
      interestGene = ranktable$UniprotID,
      referenceGene = ddd$UniprotID,
      interestGeneType = "uniprotswissprot",
      referenceGeneType = "uniprotswissprot",
      outputDirectory = this.contrast,
      isOutput = TRUE,
      perNum = nperm,
      projectName = "GSEA_proj"
    )

}


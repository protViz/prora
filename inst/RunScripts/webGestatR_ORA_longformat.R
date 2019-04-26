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

con_col <- "lhs"

contrs <- ddd %>%
  distinct(!!sym(con_col)) %>%
  pull()


# Functions ---------------------------------------------------------------

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

runWebGestaltORAlong <- function(contrast,
                                 organism = "hsapiens",
                                 ID_col = "UniprotID",
                                 target = "geneontology_Biological_Process",
                                 map_col = "GO",
                                 threshold = 0.5,
                                 direction = "greater.than",
                                 nperm = 10,
                                 contrast_col = con_col,
                                 fc_col = "estimate",
                                 outdir = "WebGestalt_ORA") {
  fpath <- make.names(contrast)

  if(!dir.exists(outdir)){
    dir.create(outdir)
  }

  dat <- ddd %>%
    dplyr::filter(!!sym(contrast_col) == contrast) %>%
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
      outputDirectory = outdir,
      isOutput = TRUE,
      perNum = nperm,
      projectName = fpath
    )
}


# Run ---------------------------------------------------------------------

sapply(contrs, runWebGestaltORAlong)

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
outdir <- "GSEA"

contrs <- ddd %>%
  distinct(!!sym(contrast_col)) %>%
  pull()


# Run ---------------------------------------------------------------------

sapply(contrs[1],
       runGSEAlong,
       organism = organism,
       ID_col = ID_col,
       target = target,
       map_col = map_col,
       nperm = nperm,
       fc_col = fc_col,
       contrast_col = contrast_col,
       outdir = outdir)

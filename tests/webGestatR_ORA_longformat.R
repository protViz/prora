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
organism <- "hsapiens"
ID_col <- "UniprotID"
target = "geneontology_Biological_Process"
map_col = "GO"
threshold = 0.5
greater = TRUE
nperm = 10
fc_col = "estimate"
outdir = "WebGestalt_ORA"

contrs <- ddd %>%
  distinct(!!sym(con_col)) %>%
  pull()


# Run ---------------------------------------------------------------------

sapply(
  contrs[1],
  runWebGestaltORAlong,
  contrast_col = con_col,
  organism = organism,
  ID_col = ID_col,
  target = target,
  map_col = map_col,
  threshold = threshold,
  greater = greater,
  nperm = nperm,
  fc_col = fc_col,
  outdir = outdir
)


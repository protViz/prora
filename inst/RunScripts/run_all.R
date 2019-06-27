rm(list = ls())
library(WebGestaltR)
library(tidyverse)
library(org.Hs.eg.db)
library(conflicted)
library(sigora)
library(GO.db)
library(slam)
library(fgczgseaora)

# Setup -------------------------------------------------------------------

# fpath <- "inst/example_data/Contrasts_SignificanceValues_f_Cells_Treatment.csv"
# dd <- read_csv(fpath)

dd <- contrast_data_example
colnames(dd) <- make.names(colnames(dd))
ddd <- getUniprotFromFastaHeader(dd)

organism <- "hsapiens"
ID_col <- "UniprotID"
target_WebGestaltR <- "geneontology_Biological_Process"
target_sigora <- "GO"
map_col <- "GO"
nperm <- 10
fc_col <- "estimate"
contrast_col <- "lhs"
threshold = 0.5
greater = TRUE
con_col <- "lhs"
outdir_GSEA <- "GSEA"
outdir_sigora <- "sigORA"
outdir_WebGestaltR = "WebGestalt_ORA"

contrs <- ddd %>%
  distinct(!!sym(contrast_col)) %>%
  pull()


# Run ---------------------------------------------------------------------

sapply(
  contrs[1],
  runGSEAlong,
  organism = organism,
  ID_col = ID_col,
  target = target_WebGestaltR,
  map_col = map_col,
  nperm = nperm,
  fc_col = fc_col,
  contrast_col = contrast_col,
  outdir = outdir_GSEA
)


# Run ---------------------------------------------------------------------

sapply(
  contrs[1],
  runSIGORAlong,
  target = target_sigora,
  fc_col = fc_col,
  greater = greater,
  fc_threshold = threshold,
  outdir = outdir_sigora
)


# Run ---------------------------------------------------------------------

sapply(
  contrs[1],
  runWebGestaltORAlong,
  contrast_col = con_col,
  organism = organism,
  ID_col = ID_col,
  target = target_WebGestaltR,
  map_col = map_col,
  threshold = threshold,
  greater = greater,
  nperm = nperm,
  fc_col = fc_col,
  outdir = outdir_WebGestaltR
)

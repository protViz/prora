# Run Script for long format

rm(list = ls())

library(org.Hs.eg.db)
library(sigora)
library(tidyverse)
library(GO.db)
library(slam)
library(fgczgseaora)


# Setup -------------------------------------------------------------------

# fpath <- "inst/example_data/Contrasts_SignificanceValues_f_Cells_Treatment.csv"
# dd <- read_csv(fpath)

dd <- contrast_data_example
colnames(dd) <- make.names(colnames(dd))

ddd <- getUniprotFromFastaHeader(dd)

con_col <- "lhs"
greater <- TRUE # TRUE: Greater than threshold, FALSE: Smaller than threshold
target <- "GO"
fc_col <- "estimate"
fc_threshold <- 0.5
outdir <- "sigORA"

contrs <- ddd %>%
  distinct(!!sym(con_col)) %>%
  pull()


# Run ---------------------------------------------------------------------

sapply(contrs[1],
       runSIGORAlong,
       target = target,
       fc_col = fc_col,
       greater = greater,
       fc_threshold = fc_threshold,
       outdir = outdir)


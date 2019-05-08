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

con_col = "lhs"
thdirection <- TRUE # TRUE: Greater than threshold, FALSE: Smaller than threshold

contrs <- ddd %>%
  distinct(!!sym(con_col)) %>%
  pull()


# Run ---------------------------------------------------------------------

sapply(contrs[1], runSIGORAlong)


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


# Run ---------------------------------------------------------------------

sapply(contrs[1], runWebGestaltORAlong)

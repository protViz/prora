rm(list = ls())
library(WebGestaltR)
library(tidyverse)
library(org.Hs.eg.db)
library(sigora)
library(GO.db)
library(slam)
library(fgczgseaora)


# Files -------------------------------------------------------------------

dd <- fgczgseaora::exampleContrastData
dfiles <- "example_data.txt"


# Parameters --------------------------------------------------------------

organism <- "hsapiens"
ID_col <- "UniprotID"
fc_col <- "estimate"

target_SIGORA <- c("GO", "KEGG", "reactome")
target_SIGORA <- target_SIGORA[1]
fc_threshold <- 0.5
greater <- TRUE


# Run ---------------------------------------------------------------------
fpath_se <- tools::file_path_sans_ext(basename(dfiles))
odir <- make.names(fpath_se)

if (!dir.exists(odir)) {
  if (!dir.create(odir)) {
    stop("\n can't create odir", odir, "\n")
  }
}

colnames(dd) <- make.names(colnames(dd))
filtered_dd <- getUniprotFromFastaHeader(dd) %>%
  filter(!is.na(UniprotID))



if (organism == "hsapiens"){

  sapply(target_SIGORA, function(target_SIGORA) {
    fgczgseaora:::.runSIGORA(
      data = filtered_dd,
      target = target_SIGORA,
      fc_col = fc_col,
      fc_threshold = fc_threshold,
      greater = greater,
      outdir = file.path(odir, "sigORA")
    )
  })
}

# Run Script for long format

rm(list = ls())

library(org.Hs.eg.db)
library(sigora)
library(tidyverse)
library(GO.db)
library(slam)
library(fgczgseaora)

fpath <- "inst/example_data/Contrasts_SignificanceValues_f_Cells_Treatment.csv"

dd <- read_csv(fpath)
colnames(dd) <- make.names(colnames(dd))

ddd <- getUniprotFromFastaHeader(dd)

trgt <- "GO"
contrast_col <- "lhs"
fc_col <- "estimate"
fc_threshold <- 0.5

contrs <- ddd %>%
  distinct(!!sym(contrast_col)) %>%
  pull()

for (this.contrast in contrs) {

  fpath <- make.names(this.contrast)

  if (!(dir.exists(fpath))) {
    dir.create(fpath)
  }

  dat <- ddd %>%
    dplyr::filter(!!sym(contrast_col) == this.contrast)

  GPStab <-
    makeGPS_wrappR(dat$UniprotID, target = trgt, dev = TRUE)

  myGPSrepo <-
    makeGPS_wrappR(dat$UniprotID, target = trgt)

  sigora_res <-
    sigoraWrappR(
      fc_threshold = fc_threshold,
      fc_col = fc_col,
      df = dat,
      GPSrepos = myGPSrepo,
      db = trgt,
      greater_than = TRUE
    )

  p1 <- try(sigora_heatmap(sigora_example, GPStab))

  rmarkdownPath <- file.path(fpath, "sigora.Rmd")

  bibpath <- file.path(fpath, "bibliography.bib")

  file.copy(
    file.path(find.package("fgczgseaora"), "rmarkdown_reports/sigora.Rmd"),
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
    params = list(results = sigora_res, plot1 = p1, GPStable = GPStab),
    clean = TRUE
  )

}


# Run Script for long format

rm(list = ls())

library(org.Hs.eg.db)
library(sigora)
library(tidyverse)
library(GO.db)
library(slam)
library(fgczgseaora)


# Setup -------------------------------------------------------------------

fpath <- "inst/example_data/Contrasts_SignificanceValues_f_Cells_Treatment.csv"

dd <- read_csv(fpath)
colnames(dd) <- make.names(colnames(dd))

ddd <- getUniprotFromFastaHeader(dd)

con_col = "lhs"

contrs <- ddd %>%
  distinct(!!sym(con_col)) %>%
  pull()


# Function ----------------------------------------------------------------

runSIGORAlong <-
  function(contrast,
           trgt = "GO",
           fc_col = "estimate",
           fc_threshold = 0.5,
           outdir = "sigORA",
           contrast_col = con_col) {

    fpath <- make.names(contrast)

    if (!dir.exists(outdir)) {
      dir.create(outdir)
    }

    if (!dir.exists(file.path(outdir, fpath))) {
      dir.create(file.path(outdir, fpath))
    }

    dat <- ddd %>%
      dplyr::filter(!!sym(contrast_col) == contrast)

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

    rmarkdownPath <- file.path(outdir, fpath, "sigora.Rmd")

    bibpath <- file.path(outdir, fpath, "bibliography.bib")

    file.copy(
      file.path(
        find.package("fgczgseaora"),
        "rmarkdown_reports/sigora.Rmd"
      ),
      rmarkdownPath,
      overwrite = TRUE
    )

    file.copy(
      file.path(
        find.package("fgczgseaora"),
        "rmarkdown_reports/bibliography.bib"
      ),
      bibpath,
      overwrite = TRUE
    )

    rmarkdown::render(
      rmarkdownPath,
      bookdown::html_document2(number_sections = FALSE),
      params = list(
        results = sigora_res,
        plot1 = p1,
        GPStable = GPStab
      ),
      clean = TRUE
    )
  }


# Run ---------------------------------------------------------------------

sapply(contrs, runSIGORAlong)


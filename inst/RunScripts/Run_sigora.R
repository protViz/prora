# Run script for sigora

rm(list = ls())

library(org.Hs.eg.db)
library(sigora)
library(tidyverse)
library(GO.db)
library(slam)
library(fgczgseaora)

fpath <-
  "inst/example_data/Ex1_interactions.csv"

dd <- read_csv(fpath)
colnames(dd) <- make.names(colnames(dd))

ddd <- getSymbolFromFasta(dd)

trgt <- "GO"

GPStab <-
  makeGPS_wrappR(ddd$Symbol, target = trgt, dev = TRUE)

myGPSrepo <-
  makeGPS_wrappR(ddd$Symbol, target = trgt) # Produce GPS repository from background

sigora_example <-
  sigoraWrappR(
    input.file = fpath,
    fc_threshold = 0.3,
    # fold change threshold of which proteins to consider "differentially expressed"
    fc_col = "estimate.Age.class..Old...Young",
    df = ddd,
    GPSrepos = myGPSrepo,
    db = trgt
  )

p1 <- sigora_heatmap(sigora_example, GPStab)


#usethis::use_data(sigora_example, overwrite = TRUE)
#usethis::use_data(p1, overwrite = TRUE)
#usethis::use_data(GPStab, overwrite = TRUE)

rmarkdown::render(
  "inst/rmarkdown_reports/sigora.Rmd",
  bookdown::html_document2(number_sections = FALSE),
  params = list(results = sigora_example, plot1 = p1, GPStable = GPStab),
  clean = TRUE
)

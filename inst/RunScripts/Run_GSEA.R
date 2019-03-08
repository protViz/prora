# Run script for GSEA using WebGestaltR

rm(list = ls())
library(WebGestaltR)
library(tidyverse)
library(fgczgseaora)
library(conflicted)

fpath <-
  "inst/example_data/Ex1_interactions.csv"

dd <- read_csv(fpath)
colnames(dd) <- make.names(colnames(dd))
ddd <- getSymbolFromSwissprotID(dd)

outputDir <- "GSEA_output"
dir.create(outputDir)

organism <- "hsapiens"
ID_col <- "Symbol"
fc_col <- "estimate.Age.class..Old...Young"
target <- "geneontology_Biological_Process"

input <- ddd[, c(ID_col, fc_col)]

GSEA_res <-
  WebGestaltR(
    enrichMethod = "GSEA",
    # does permutation test with 1000 permutations per default, might take a while
    organism = organism,
    enrichDatabase = target,
    interestGene = input,
    interestGeneType = "genesymbol",
    # Or "uniprotswissprot" in some cases
    outputDirectory = "GSEA_output/",
    isOutput = TRUE,
    perNum = 1000,
    projectName = "GSEA_proj"
  )

GSEA <- list(
  organism = organism,
  target = target,
  outputDir = outputDir,
  input_data = input,
  GSEA_res = GSEA_res
)

#usethis::use_data(GSEA, overwrite = TRUE)

rmarkdown::render(
  "inst/rmarkdown_reports/GSEA.Rmd",
  bookdown::html_document2(number_sections = FALSE),
  params = list(GSEA = GSEA),
  clean = TRUE
)

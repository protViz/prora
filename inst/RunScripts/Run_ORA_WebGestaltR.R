rm(list = ls())
library(WebGestaltR)
library(tidyverse)
library(fgczgseaora)
library(org.Hs.eg.db)
library(conflicted)

fpath <- "inst/example_data/Ex1_interactions.csv"

dd <- read_csv(fpath)
colnames(dd) <- make.names(colnames(dd))
dd <- getUniprotFromFastaHeader(dd)

organism <- "hsapiens"
ID_col <- "UniprotID"
target <- "geneontology_Biological_Process"
map_col <- "GO"
threshold <- log2(1.5)
direction <- "greater.than"
nperm <- 10

columns <- grep("estimate", colnames(dd), value = TRUE)

fc_col <- columns[2]

if (!dir.exists(fc_col)) {
  dir.create(fc_col)
}

apply_threshold <-
  function(df,
           th,
           fc_col,
           ID_col,
           alt = c("greater.than", "less.than")) {
    alt <- match.arg(alt)
    if (alt == "greater.than") {
      df %>%
        dplyr::filter(!!sym(fc_col) > th) %>%
        dplyr::select(!!sym(ID_col), !!sym(fc_col)) -> out
    } else if (alt == "less.than") {
      df %>%
        dplyr::filter(!!sym(fc_col) <= th) %>%
        dplyr::select(!!sym(ID_col), !!sym(fc_col)) -> out
    }
    return(out)
  }

ranktable <-
  apply_threshold(
    df = dd,
    th = threshold,
    fc_col = fc_col,
    ID_col = ID_col,
    alt = direction
  )

ORA_res <-
  WebGestaltR(
    enrichMethod = "ORA",
    organism = organism,
    enrichDatabase = target,
    interestGene = ranktable$UniprotID,
    referenceGene = dd$UniprotID,
    interestGeneType = "uniprotswissprot",
    referenceGeneType = "uniprotswissprot",
    outputDirectory = fc_col,
    isOutput = TRUE,
    perNum = nperm,
    projectName = "GSEA_proj"
  )

rm(list = ls())
library(WebGestaltR)
library(tidyverse)
library(fgczgseaora)
library(org.Hs.eg.db)
library(conflicted)

fpath <- "inst/example_data/Ex1_interactions.csv"

dd <- read_csv(fpath)
colnames(dd) <- make.names(colnames(dd))
dd <- getUniprotFromFastaHeader(dd) %>%
  dplyr::filter(!is.na(UniprotID))

organism <- "hsapiens"
ID_col <- "UniprotID"
target <- "geneontology_Biological_Process"
map_col <- "GO"
nperm <- 10

columns <- grep("estimate",colnames(dd), value=TRUE)

fc_col <- columns[2]

if(!dir.exists(fc_col)){
  dir.create(fc_col)
}

ranktable <-
  dplyr::select(dd, ID = !!sym(ID_col), Score = !!sym(fc_col))

GSEA_res <-
  WebGestaltR(
    enrichMethod = "GSEA",
    organism = organism,
    enrichDatabase = target,
    interestGene = ranktable,
    interestGeneType = "uniprotswissprot",
    outputDirectory = fc_col,
    isOutput = TRUE,
    perNum = nperm,
    projectName = "GSEA_proj"
  )

mappingTable <-
  read_delim(file.path(fc_col, "Project_GSEA_proj/interestingID_mappingTable_GSEA_proj.txt"),
             delim = "\t")
mappingTable %>% mutate(entrezgene = as.character(entrezgene)) -> mappingTable
GSEA_res_sep <- GSEA_res %>% separate_rows(leadingEdgeId, sep=";")
merged_data <- inner_join(mappingTable, GSEA_res_sep,
                          by=c("entrezgene" = "leadingEdgeId"))

GSEA <- list(
  organism = organism,
  target = target,
  input_data = ranktable,
  output_dir = fc_col,
  merged_data = merged_data,
  nperm = nperm
)

rmarkdownPath <- file.path(fc_col, "GSEA.Rmd")
bibpath <- file.path(fc_col, "bibliography.bib")
file.copy(
  file.path(find.package("fgczgseaora"), "rmarkdown_reports/GSEA.Rmd"),
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
  params = list(GSEA = GSEA),
  clean = TRUE
)

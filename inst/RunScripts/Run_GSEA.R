# Run script for GSEA using WebGestaltR

rm(list = ls())
library(WebGestaltR)
library(tidyverse)
library(fgczgseaora)
library(org.Hs.eg.db)
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
# fc_col <- "estimate.Age.class..Old...Young"
fc_col <- "estimate.Age.classOld..PathologyHF.iDCM...PathologyNone"
target <- "geneontology_Biological_Process"
map_col <- "GO"

input <- ddd[ , c(ID_col, fc_col)]
colnames(input) <- c("ID", "Score")

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

mappingTable <- read_delim("GSEA_output/Project_GSEA_proj/interestingID_mappingTable_GSEA_proj.txt", delim = "\t")
mappingTable$entrezgene <- as.character(mappingTable$entrezgene)

merged_data <- AnnotationDbi::mapIds(org.Hs.eg.db,
                             as.character(mappingTable$entrezgene),
                             column = map_col,
                             keytype = "ENTREZID",
                             multiVals = "CharacterList") %>%
  unlist() %>% data.frame(entrezgene = names(.), geneSet = .) %>%
  dplyr::inner_join(., GSEA_res) %>%
  dplyr::inner_join(., mappingTable)

toplot <- merged_data %>%
  dplyr::select(geneSet, geneSymbol) %>%
  dplyr::mutate(IDD = 1:nrow(.)) %>%
  tidyr::spread(geneSet, geneSymbol) %>%
  dplyr::select(-IDD) %>%
  as.list() %>%
  lapply(na.omit)


GSEA <- list(
  organism = organism,
  target = target,
  outputDir = outputDir,
  input_data = input,
  GSEA_res = GSEA_res,
  merged_data = merged_data,
  upsetData = toplot
)

# usethis::use_data(GSEA, overwrite = TRUE)

rmarkdown::render(
  "inst/rmarkdown_reports/GSEA.Rmd",
  bookdown::html_document2(number_sections = FALSE),
  params = list(GSEA = GSEA),
  clean = TRUE
)

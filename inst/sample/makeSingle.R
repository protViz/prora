rm(list = ls())
library(fgsea)
library(msigdbr)
library(tidyverse)
library(fgcz.gsea.ora)

#data <- readxl::read_xlsx("d:/Dropbox/DataAnalysis/p3433_o7341_20200917/Contrasts_GLM.xlsx")
if (TRUE) {

  data <- readxl::read_xlsx("p3433_results/ImputationResults_p3433_model_1_All/contrastsWithImputation_all.xlsx")
  unique(data$contrast)
  data2 <- fgcz.gsea.ora::get_UniprotID_from_fasta_header(data)
  prefix <- "Imputation_"
  data3 <- fgcz.gsea.ora::map_ids_uniprot(data2)
  ranklist <- fgsea_rank_contrasts(data3, ids = "P_ENTREZGENEID", score = "estimate_median"  )

}else{

  data <- readxl::read_xlsx("p3433_results/LinearModelResult_p3433_model_1_noKO/Contrasts_Model_B.xlsx")
  data2 <- data
  prefix <- "LinearModel_"
  data3 <- fgcz.gsea.ora::map_ids_uniprot(data2)
  ranklist <- fgsea_rank_contrasts(data3, ids = "P_ENTREZGENEID", score = "statistic"  )
}


msigdbr::msigdbr_species()
species <- "Homo sapiens"

C5 <- msigdbr_collections() %>% filter(gs_cat == "C5")
fgseaGSlist <- fgsea_msigdb_collections(C5, species = "Homo sapiens")

i <- 1
print(names(fgseaGSlist)[i])
fgseaGS <- fgseaGSlist[[i]]

fgseaRes <- run_fgsea_for_allContrasts(ranklist,fgseaGS, nperm = 10000  )
allres <- dplyr::bind_rows(fgseaRes)

writexl::write_xlsx(allres, path = paste0(prefix, names(fgseaGSlist)[i], ".xlsx"))

threshold <- 0.1

for (contrast in 1:length(fgseaRes)) {
  print(contrast)
  fgseaResult <- fgseaRes[[contrast]]
  relevantResult <- fgseaResult %>% dplyr::filter(padj < threshold)
  geneSet <- fgseaGS
  rankList <- ranklist[[contrast]]
  dim(relevantResult)
  GSEAResults <- list(
    threshold = threshold,
    gsName = names(fgseaGSlist)[i],
    fgseaRes = fgseaResult,
    geneSet = geneSet,
    rankList = rankList,
    relevantResult = relevantResult
  )

  outfile = paste0(prefix , make.names(names(ranklist)[contrast]), names(fgseaGSlist)[i], ".html")
  rmarkdown::render("VisualizeSingle.Rmd", params = list(GSEAResults = GSEAResults), output_file = outfile)
}

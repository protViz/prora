rm(list = ls())
library(fgsea)
library(msigdbr)
library(tidyverse)
library(fgcz.gsea.ora)

useLog2FC <- TRUE
inputData <- "2GrpAppOutput/MQ_2grp_report_o6927_iMPC_vs_Ctrl.txt"
#species <- "Homo sapiens"
species <- "Mus musculus"

outname <- gsub("\\.txt","",basename(inputData))


#####################
#########################



padj_threshold <- 0.1
prefix <- "FGSEA_"
msigdbr::msigdbr_species()

summaries <- list()
data <- readr::read_tsv(inputData)
summaries$nrow_all <- nrow(data)
# remove reverse sequences and contaminants
data <- data %>% filter(!(grepl("^REV__", TopProteinName ) | grepl("^zz\\|", TopProteinName)))

# Extract uniprotIDS
data2 <- fgcz.gsea.ora::get_UniprotID_from_fasta_header(data,
                                                        idcolumn = "TopProteinName")

# map to entriz id's
data3 <- fgcz.gsea.ora::map_ids_uniprot(data2)

ranklist <- data3 %>% dplyr::select(dplyr::all_of(c("P_ENTREZGENEID","pseudo.log2FC", "pseudo.P.Value" )))
ranklist <- na.omit(ranklist)

if (useLog2FC) {
  ranklist <- dplyr::select(ranklist, P_ENTREZGENEID, score = pseudo.log2FC)
} else {
  ranklist$score <- sign(ranklist$pseudo.log2FC) * (1 - ranklist$pseudo.P.Value)
  ranklist <- dplyr::select(ranklist, dplyr::all_of(c("P_ENTREZGENEID","score")))
}

ranklist <- ranklist %>% group_by(!!sym("P_ENTREZGENEID")) %>% summarize(score = mean(score),.groups = "drop" )
summaries$nrow_mapped <- nrow(ranklist)
readr::write_tsv(ranklist, path =  paste0(prefix, outname, ".rnk"), col_names = FALSE)

rankarray <- ranklist$score
names(rankarray) <- ranklist$P_ENTREZGENEID

#### Run GSEA analysis


C5 <- msigdbr_collections() %>% filter(gs_cat == "C5")
C5 <- C5[1:3,]
fgseaGSlist <- fgsea_msigdb_collections(C5, species = species)
names(C5)


#fgsea(pathways =  fgseaGSlist[[1]], rankarray)


fgseaRes <- run_fgsea_for_allGeneSets(rankarray, fgseaGSlist, nperm = 10000  )
allres <- dplyr::bind_rows(fgseaRes)

rankList <- rankarray


for (iGS in 1:length(fgseaRes)) {
  iGS <- 1
  fgseaResult <- fgseaRes[[iGS]]
  geneSet <- fgseaGSlist[[iGS]]

  relevantResult <- fgseaResult %>% dplyr::filter(padj < padj_threshold)


  collapsedPathways <- collapsePathways(relevantResult,
                                        geneSet, rankList)
  relevantResult <- dplyr::select(relevantResult, -ES)

  mainPathways <- fgseaResult %>% dplyr::filter(pathway %in% collapsedPathways$mainPathways) %>%
    dplyr::select( -ES)

  GSEAResults <- list(
    threshold = padj_threshold,
    gsName = names(fgseaGSlist)[iGS],
    fgseaRes = fgseaResult,
    geneSet = geneSet,
    rankList = rankList,
    relevantResult = relevantResult,
    mainPathways = mainPathways,
    map_summaries = summaries,
    score = if (useLog2FC) {"log2FC"} else {"scaled p-value"}
  )

  outfile = paste0(outname , names(fgseaGSlist)[iGS])
  html_out <- paste0(prefix, outfile, ".html")
  rmarkdown::render("VisualizeSingle.Rmd",
                    params = list(GSEAResults = GSEAResults),
                    output_file = html_out)
  writexl::write_xlsx(mainPathways, path = paste0(prefix,  outfile, "_MainPathways" , ".xlsx"))
  writexl::write_xlsx(fgseaResult, path = paste0(prefix, outfile, "_AllPathways", ".xlsx" ))
}



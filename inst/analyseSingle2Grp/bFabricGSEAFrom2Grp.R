#! /usr/lib/Rscript

library(fgsea)
library(msigdbr)
library(tidyverse)
library(prora)
library(yaml)
library(writexl)

YAML = TRUE
useLog2FC <- TRUE

if (YAML) {
  args = commandArgs(trailingOnly = TRUE)
  if (length(args) > 0) {
    yamlfile <- args[1]
  }else{
    stop("script needs one argument the bfabripy.yaml file.")
    yamlfile <- "WU256806.yaml"
  }
  parameters <- yaml::read_yaml(yamlfile)

  inputData <- basename(parameters$application$input$`MaxQuant - Two Group Analysis Report`)
  FDR_threshold <- as.numeric(parameters$application$parameters$FDRthreshold)
  species <- parameters$application$parameters$Species
  outname <- gsub("\\.txt","",basename(inputData))

} else {

  inputData <- "2GrpAppOutput/MQ_2grp_report_o6927_iMPC_vs_Ctrl.txt"
  inputData <- "2GrpAppOutput/MQ_2grp_report_o6927_MyoD_D10_vs_Ctrl.txt"
  inputData <- "2GrpAppOutput/MQ_2grp_report_o6927_MyoD_FRC_D10_vs_Ctrl.txt"
  #inputData <- "2GrpAppOutput/MQ-report-bfabric-FASPvsUrea_jg_n_EAR.zip.txt"

  species <- "Homo sapiens"
  #species <- "Mus musculus"
  FDR_threshold <- 0.1 # (0.05, 0.1, 0.25)
  workUnit <- "WU300000"
  outname <- gsub("\\.txt","",basename(inputData))
}


outdir <- "out_dir"

prora::copy_bfabric_2grp()

#if (dir.exists(outdir)) {
#  unlink(outdir,recursive = TRUE)
#}

prefix <- "FGSEA_"

#####################

stopifnot(species %in% msigdbr::msigdbr_species()$species_name)
dir.create(outdir)


summaries <- list()
data <- readr::read_tsv(inputData)

# remove reverse sequences and contaminants
data <- data %>% filter(!(grepl("^REV__", TopProteinName ) | grepl("^zz\\|", TopProteinName)))
summaries$nrow_all <- nrow(data)

# Extract uniprotIDS
data2 <- prora::get_UniprotID_from_fasta_header(data,
                                                        idcolumn = "TopProteinName")
# map to entriz id's
data3 <- prora::map_ids_uniprot(data2)
mappingtable <- data3 %>% dplyr::select(all_of(c("UniprotID", "P_ENTREZGENEID", "TopProteinName")))

writexl::write_xlsx(mappingtable,
                    path = file.path(outdir,paste0(prefix, "Mapping_" , outname,".xlsx")))

### prepare ranklist
ranklist <- data3 %>% dplyr::select(dplyr::all_of(c("P_ENTREZGENEID","pseudo.log2FC", "pseudo.P.Value" )))
ranklist <- na.omit(ranklist)

if (useLog2FC) {
  ranklist <- dplyr::select(ranklist, .data$P_ENTREZGENEID, score = .data$pseudo.log2FC)
} else {
  ranklist$score <- sign(ranklist$pseudo.log2FC) * (1 - ranklist$pseudo.P.Value)
  ranklist <- dplyr::select(ranklist, dplyr::all_of(c("P_ENTREZGENEID","score")))
}

ranklist <- ranklist %>%
  dplyr::group_by(.data$P_ENTREZGENEID) %>%
  dplyr::summarize(score = mean(.data$score),.groups = "drop" )

if (nrow(ranklist) == 0) {
  ErrorMessage <- "No Id's were mapped."
  rmarkdown::render("ErrorMessage.Rmd",
                    params = list(ErrorMessage = ErrorMessage,
                                  protIDs = sample(data$TopProteinName, size = 10)))
  file.copy("ErrorMessage.html", file.path(outdir, "ErrorMessage.html") , overwrite = TRUE)

  write(ErrorMessage,
        file = stderr())
  stop(ErrorMessage)
}

summaries$nrow_mapped <- nrow(ranklist)

readr::write_tsv(ranklist,
                 file =  file.path(outdir, paste0(prefix, outname, ".rnk")),
                 col_names = FALSE)

rankarray <- ranklist$score
names(rankarray) <- ranklist$P_ENTREZGENEID

#### Run GSEA analysis

C5 <- msigdbr::msigdbr_collections() %>%
  filter(gs_cat == "C5" | gs_cat == "H")

C5 <- C5[c(1:3,5), ]

fgseaGSlist <- fgsea_msigdb_collections(C5, species = species)

#fgsea(pathways =  fgseaGSlist[[1]], rankarray)
#undebug(run_fgsea_for_allGeneSets)
fgseaRes <- run_fgsea_for_allGeneSets(rankarray, fgseaGSlist, nperm = 10000  )

allres <- dplyr::bind_rows(fgseaRes)

if (nrow(allres) == 0) {
  # write HTML
  ErrorMessage <- paste0("No results for all GeneSets, please do check if species ",
                         species, " is correct.\n")
  rmarkdown::render("ErrorMessage.Rmd",
                    params = list(ErrorMessage = ErrorMessage,
                                  protIDs = sample(data$TopProteinName, size = 10)))


  file.copy("ErrorMessage.html", file.path(outdir, "ErrorMessage.html") , overwrite = TRUE)

  write(ErrorMessage,
        file = stderr())
  stop(ErrorMessage)
} else{
  writexl::write_xlsx(allres,
                      path = file.path(outdir, paste0(prefix,"All_results",outname , ".xlsx")))
}


select_relevant_results <- function(fgseaResult,
                                    geneSet,
                                    gsName,
                                    rankList,
                                    FDR_threshold){
  relevantResult <- fgseaResult %>% dplyr::filter(.data$padj < FDR_threshold)
  collapsedPathways <- fgsea::collapsePathways(relevantResult,
                                               geneSet, rankList)
  relevantResult <- dplyr::select(relevantResult, -ES)
  mainPathways <- fgseaResult %>% dplyr::filter(.data$pathway %in% collapsedPathways$mainPathways) %>%
    dplyr::select( -ES)

  GSEAResults <- list(
    threshold = FDR_threshold,
    gsName = gsName,
    fgseaRes = fgseaResult,
    geneSet = geneSet,
    rankList = rankList,
    relevantResult = relevantResult,
    mainPathways = mainPathways,
    map_summaries = summaries,
    score = if (useLog2FC) {"log2FC"} else {"scaled p-value"}
  )
  return(GSEAResults)
}


for (iGS in 1:length(fgseaRes)) {
  #iGS <- 1
  fgseaResult <- fgseaRes[[iGS]]
  geneSet <- fgseaGSlist[[iGS]]
  gsName = names(fgseaGSlist)[iGS]

  if (!is.null(fgseaResult)) {
    #undebug(select_relevant_results)
    GSEAResults <- select_relevant_results(fgseaResult,
                                           geneSet,
                                           gsName,
                                           rankarray,
                                           FDR_threshold)




    outfile = paste0(outname ,gsName)
    html_out <- paste0(prefix, outfile, ".html")
    rmarkdown::render("VisualizeSingle.Rmd",
                      params = list(GSEAResults = GSEAResults))
    file.copy("VisualizeSingle.html", file.path(outdir, html_out) , overwrite = TRUE)

    writexl::write_xlsx(GSEAResults$mainPathways,
                        path = file.path(outdir,
                                         paste0(prefix,  outfile, "_MainPathways" , ".xlsx")))
    writexl::write_xlsx(fgseaResult,
                        path = file.path(outdir,
                                         paste0(prefix, outfile, "_AllPathways", ".xlsx" )))

  } else {
    write(paste0("No results for GS : ",gsName, "\n"), file = stderr())
  }

}



# Author Witold Wolski
# RscriptExe <- "c:/Program Files/R/R-4.1.0/bin/Rscript.exe"
# arguments <- c(rscriptlocation, path, "human", parameters$outfolder , parameters$clustalg , workunitid, projectid, parameters$peptide, parameters$JK)
# res <- system2(RscriptExe, args = arguments)

library(tidyverse)
library(prolfqua)
#library(dendextend)
library(dynamicTreeCut)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(prora)

prora::copy_bfabric_ClusterProfiler()

#args = commandArgs(trailingOnly = TRUE)
args <- list()
args$l1 <- 1
parameter <- list()

parameter$nrCluster <- 10
parameter$row_scale <- TRUE
parameter$row_center <- TRUE
parameter$transform = "log2scaled" # "log2"
parameter$peptide  = FALSE
parameter$JK = TRUE


map_species <- function(thisname){
  if (thisname == "Homo sapiens") {return("human")}
  if (thisname == "Mus musculus") {return("mouse")}
}
if ( length(args) == 0 ) {
  # parameters
  datadir <- file.path(find.package("prolfquaData") , "quantdata")
  parameter$inputMQfile <-  file.path(datadir, "MAXQuant_ComboCourse_p2691_March_2018_WU183012.zip")
  parameter$organism <- "yeast" # "human", "mouse"
  parameter$outpath = "dummy"
  parameter$clustering <- "DPA" #
  #parameter$clustering <- "hclustdeepsplit"
  parameter$projectID <- 3000
  parameter$workunitID <- 233333
} else if ( length(args) == 1) {
  parameters <- yaml::read_yaml("WU268998.yaml")
  parameters <- yaml::read_yaml("WU269043.yaml")
  parameters <- yaml::read_yaml("WU269045.yaml")

  print("::::: USING YAML FILE ::::")

  parameter$inputMQfile <- basename(parameters$application$input$MaxQuant)
  parameter$organism <- map_species(parameters$application$parameters$Species)
  parameter$outpath <- "out_dir"
  parameter$pthreshold <- as.numeric(parameters$application$parameters$FDRthreshold)
  parameter$clustering <- parameters$application$parameters$ClustAlg
  parameter$workunitID <- parameters$job_configuration$workunit_id
  parameter$projectID <- gsub(".*htdocs/(p[0-9]{3,8})/bfabric.*","\\1",parameters$application$output)
  parameter$peptide <- parameters$application$parameters$MQInputFile == "peptide.txt"
  parameter$JK <- "euclidean_JK" == parameters$application$parameters$Distance
  print(parameter)
  # read yaml and extract
} else {
  print("::::USING COMMAND LINE ARGS:::")
  parameter$inputMQfile <- args[1]
  parameter$organism <- args[2]
  parameter$outpath <- args[3]
  parameter$clustering <- args[4]
  parameter$workunitID <- args[5]
  parameter$projectID <- args[6]
  parameter$peptide <- as.logical(args[7])
  parameter$JK <- as.logical(args[8])
  print(parameter)
  print("::::END OF PARAMS::::")
}


if (parameter$clustering == "DPA") {
  #remotes::install_github("mariaderrico/DPAclustR")
  library(DPAclustR)
  # Install DPA Python package from local path using reticulate
  # Prerequisite:
  #   - install the DPA package globally or in the virtualenv using
  #     the "Installation" instructions at https://github.com/mariaderrico/DPA.
  #   - specify "path_to_python_virtualenv" and "path_to_DPA_local_package" below
  library(reticulate)

  #if (dir.exists("/scratch/CLUSTERPORFILER/pythonenv1")) {
  #  reticulate::use_virtualenv("/scratch/CLUSTERPORFILER/pythonenv1")
  #}

  # Python version configuration (choose one of the three options below):
  # use_python("path_to_python_binary")
  # use_condaenv("path_to_conda_environmet")
  #use_virtualenv("path_to_python_virtualenv", required=TRUE)
  #py_install("DPA")
  use_virtualenv("/scratch/CLUSTERPROFILER/pythonenv2", required = TRUE)
  DPA <- reticulate::import("Pipeline.DPA")
}


# related to data preprocessing



# related to clustering


# RESULTS will be saved as rds file
RESULTS <- list()


if (!dir.exists(parameter$outpath)) {
  dir.create(parameter$outpath)
}

# Set up species
if (parameter$organism == "yeast") {
  orgDB <- "org.Sc.sgd.db"
  parameter$species <- "Saccharomyces cerevisiae"
} else if (parameter$organism == "human") {
  orgDB <- "org.Hs.eg.db"
  parameter$species <-  "Homo sapiens"
} else if (parameter$organism == "mouse") {
  orgDB <- "org.Mm.eg.db"
  parameter$species <- "Mus musculus"
}

stopifnot(parameter$species %in% msigdbr::msigdbr_species()$species_name)


# Read peptides or proteins
if (parameter$peptide) {
  datatmp <- prolfqua::tidyMQ_Peptides(parameter$inputMQfile)
  config <- prolfqua::create_config_MQ_peptide()
} else {
  datatmp <- prolfqua::tidyMQ_ProteinGroups(parameter$inputMQfile)
  datatmp <- datatmp %>% filter(.data$nr.peptides > 1)
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "raw.file"
  atable$hierarchy[["protein_Id"]] <- c("majority.protein.ids")
  atable$hierarchyDepth <- 1
  atable$setWorkIntensity("mq.protein.intensity")
  config <- AnalysisConfiguration$new(atable)
}

datatmp$AllSamples <- "S"
datatmp$file <- datatmp$raw.file
config$table$factors[["all"]] = "AllSamples"
config$table$factors[["file"]] = "file"



sdata <- setup_analysis(datatmp, config)
lfqd <- LFQData$new(sdata, config)
lfqd$filter_proteins_by_peptide_count()

lfqd$remove_small_intensities()
lfqd$hierarchy_counts()


# different ways to normalize data
{
  st <- lfqd$get_Stats()
  tmp <- st$stats()
  RESULTS$CV.50 <- median(tmp$CV, na.rm = TRUE)

  if (parameter$transform == "log2") {
    print("ONLY log2 transform intensities.")
    tr <- lfqd$get_Transformer()
    transformed <- tr$log2()$lfq
  } else if (parameter$transform == "log2scaled") {
    tr <- lfqd$get_Transformer()
    transformed <- tr$log2()$robscale()$lfq
  }
  if (parameter$peptide) {
    ag <- transformed$get_Aggregator()
    ag$medpolish()
    prot <- ag$lfq_agg
  } else {
    prot <- transformed
  }
  statstmp <- prot$get_Stats()$stats()
  RESULTS$SD.50 <- median(statstmp$sd, na.rm = TRUE)
}



RESULTS$prot <- prot

wide <- prot$to_wide(as.matrix = TRUE)
mdata <- wide$data

RESULTS$dataDims <- c(nrPort = nrow(mdata))


cp_filterforNA  <- function(mdata)
{
  na <- apply(mdata, 1, function(x) {
    sum(is.na(x))
  })
  nc <- ncol(mdata)
  mdata <- mdata[na < floor(0.5 * nc) - 1, ]
  return(mdata)
}



mdataf <- cp_filterforNA(mdata)
RESULTS$dataDims <- c(RESULTS$dataDims, nrPortNoNas = nrow(mdataf))

# scale matrix rows
scaledM <- t(scale(t(mdataf),center = parameter$row_center, scale = parameter$row_scale))
# drop rows with NA produced because of zero variance.
scaledM <- scaledM[!(rowSums(is.na(scaledM)) == ncol(scaledM)),]


if (parameter$clustering == "hclust") {
  resClust <- cp_clusterHClustEuclideanDist(scaledM,nrCluster = parameter$nrCluster, JK = parameter$JK)
} else if (parameter$clustering == "hclust_deepsplit") {
  resClust <- cp_clusterHClustEuclideanDistDeepslit(scaledM, JK = parameter$JK)
} else if (parameter$clustering == "DPA") {
  resClust <- cp_clusterDPAEuclideanDist(round(scaledM,10)[!duplicated(round(scaledM,10)),], metric = "precomputed", JK = parameter$JK)
}

RESULTS$scaledM <- scaledM
RESULTS$dendrogram <- resClust$dendrogram
RESULTS$nrCluster <- resClust$nrCluster


clusterAssignment <- prora::get_UniprotID_from_fasta_header(
  resClust$clusterAssignment,
  idcolumn = "protein_Id")

RESULTS$dataDims <- c(RESULTS$dataDims,  UniprotExtract = sum(!is.na(clusterAssignment$UniprotID)))


res <- map_ids_2ways(clusterAssignment, species = parameter$species)
clusterAssignment <- res$clusterAssignment
RESULTS$id.mapping.service <- res$mapping.service

# clusterProfiler ----
#clusterAssignmentF

RESULTS$dataDims <- c(RESULTS$dataDims,  ENTREZGENEID = sum(!is.na(clusterAssignment$P_ENTREZGENEID)))
RESULTS$clusterAssignment <- clusterAssignment
clusterB <- clusterAssignment %>% filter(!is.na(P_ENTREZGENEID))
# check mapping efficiency.


#clusterB$clusterID <- paste0("Cluster", clusterB$clusterID)
clusterProfilerinput <- split(clusterB$P_ENTREZGENEID, clusterB$Cluster)


resGOEnrich <- list()
mt <- "GO Biological Process"
clustProf <- cp_compClust(clusterProfilerinput,orgDB,clusterB$P_ENTREZGENEID,ont = "BP")
if (!is.character(clustProf)) {
  resGOEnrich$BP <- list(mt = mt,clustProf = clustProf)
}

mt <- "GO Molecular Function"
clustProf <- cp_compClust(clusterProfilerinput, orgDB, clusterB$P_ENTREZGENEID,ont = "MF")
if (!is.character(clustProf)) {
  resGOEnrich$MF <- list(mt = mt, clustProf = clustProf)
}

mt <- "GO Cellular Component"
clustProf <- cp_compClust(clusterProfilerinput, orgDB, clusterB$P_ENTREZGENEID,ont = "CC")
if (!is.character(clustProf)) {
  resGOEnrich$CC <- list(mt = mt, clustProf = clustProf)
}

RESULTS$resGOEnrich <- resGOEnrich

outfile <- tools::file_path_sans_ext(basename(parameter$inputMQfile))
outfile <- paste0(parameter$projectID, "_",
                  as.character(parameter$workunitID), "_",
                  as.character(outfile), "_",
                  parameter$clustering)


#saveRDS(RESULTS,
#        file = file.path(parameter$outpath, paste0(outfile,".Rds")))


## create summary.

EXCEL_RESULTS <- list()

output2 <- lapply(RESULTS$resGOEnrich,
                  function(x){res <- as.data.frame(x$clustProf); res$GS <- x$mt; res})

if (length(output2) > 0) {
  output2 <- dplyr::bind_rows(output2)

  output2 <- tibble::add_column(output2,
                                projectID = parameter$projectID,
                                workunitID = parameter$workunitID,
                                zipfile = basename(parameter$inputMQfile),
                                clustering  = parameter$clustering, .before = 1 )

  #output2 <- output2 %>% filter(p.adjust < parameter$pthreshold)
  #readr::write_tsv(output2 , file = file.path(parameter$outpath, paste0("GS_", outfile, ".tsv")))
  EXCEL_RESULTS$GS <- output2

  RESULTS$nr.of.GS <- output2 %>% filter(p.adjust < parameter$pthreshold) %>% nrow
} else {
  RESULTS$nr.of.GS <- 0
}



# output :
output1 <- data.frame(
  projectID = parameter$projectID,
  workunitID = parameter$workunitID,
  zipfile = basename(parameter$inputMQfile),
  clustering  = parameter$clustering,
  peptide = as.character(parameter$peptide),
  JK = as.character(parameter$JK),
  nr.proteins = RESULTS$dataDims["nrPort"],
  nr.proteins.NA.filtered = RESULTS$dataDims["nrPortNoNas"],
  nr.UniprotIDs = RESULTS$dataDims["UniprotExtract"],
  nr.ENTREZIDS = RESULTS$dataDims["ENTREZGENEID"],
  nr.samples = nrow(RESULTS$prot$factors()),
  nr.of.clusters = RESULTS$nrCluster,
  FDRThreshold = parameter$pthreshold,
  nr.of.GS.025 = RESULTS$nr.of.GS,
  median.CV = RESULTS$CV.50,
  median.sd = RESULTS$SD.50,
  id.mapping.service = RESULTS$id.mapping.service
)

EXCEL_RESULTS$summary <- output1




# Protein_zipfilename_clusteringname.txt
# 300 - 5000 rows
# zipfile
# clusterinng
# Protein IDS
# cluster assignments
# Protein intensities sample_Id

output3 <- RESULTS$clusterAssignment
output3$projectID = parameter$projectID
output3$workunitID = parameter$workunitID
output3$zipfile = basename(parameter$inputMQfile)
output3$clustering  = parameter$clustering
output3$peptide = as.character(parameter$peptide)
output3$JK = as.character(parameter$JK)

tmp <- RESULTS$prot$to_wide()$data
output3 <- right_join(output3, tmp, by = "protein_Id")

protfilename <- paste0("Protein_" , outfile, '.tsv')

EXCEL_RESULTS$proteinData <- output3


filermd <- paste0("tmp_profileC",paste(sample(LETTERS, 5, TRUE), collapse = ""))
file.copy("profileClusters_V2.Rmd", paste0(filermd,".Rmd"), overwrite = TRUE)

rmarkdown::render(paste0(filermd,".Rmd"),
                  params = list(resultsxx = RESULTS, parametersxx = parameter))

file.copy(paste0(filermd,".html"),
          file.path(parameter$outpath, paste0("HTML_",outfile, ".html")),overwrite = TRUE)

file.remove(paste0(filermd,".Rmd"),paste0(filermd,".html"))

writexl::write_xlsx(EXCEL_RESULTS, path = file.path(parameter$outpath, paste0("",outfile, ".xlsx")))

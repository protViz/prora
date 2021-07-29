library(tidyverse)
library(prolfqua)
#library(dendextend)
library(dynamicTreeCut)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(prora)

args = commandArgs(trailingOnly = TRUE)

parameter <- list()

if ( length(args) == 0 ) {
  # parameters
  datadir <- file.path(find.package("prolfquaData") , "quantdata")
  parameter$inputMQfile <-  file.path(datadir, "MAXQuant_ComboCourse_p2691_March_2018_WU183012.zip")
  parameter$organism <- "yeast" # "human", "mouse"
  parameter$outpath = "dummy"
  parameter$clustering <- "DPA"
  #parameter$clustering <- "hclustdeepsplit"
  parameter$projectID <- 3000
  parameter$workunitID <- 233333
} else if ( length(args) == 1) {
  # read yaml and extract
} else {
  print("::::USING COMMAND LINE ARGS:::")
  parameter$inputMQfile <- args[1]
  parameter$organism <- args[2]
  parameter$outpath <- args[3]
  parameter$clustering <- args[4]
  parameter$workunitID <- args[5]
  parameter$projectID <- args[6]
  print(parameter)
}


if (parameter$clustering == "DPA") {
  #remotes::install_github("mariaderrico/DPAclustR")
  library(DPAclustR)
  # Install DPA Python package from local path using reticulate
  # Prerequisite:
  #   - install the DPA package globally or in the virtualenv using
  #     the "Installation" instructions at https://github.com/mariaderrico/DPA.
  #   - specify "path_to_python_virtualenv" and "path_to_DPA_local_package" below

  #if (dir.exists("/scratch/CLUSTERPORFILER/pythonenv1")) {
  #  reticulate::use_virtualenv("/scratch/CLUSTERPORFILER/pythonenv1")
  #}

  # Python version configuration (choose one of the three options below):
  # use_python("path_to_python_binary")
  # use_condaenv("path_to_conda_environmet")
  #use_virtualenv("path_to_python_virtualenv", required=TRUE)
  #py_install("DPA")
  DPA <- reticulate::import("Pipeline.DPA")
}


# related to data preprocessing
parameter$peptide  = FALSE
parameter$transform = "log2scaled" # "log2"

parameter$row_center <- TRUE
parameter$row_scale <- TRUE

# related to clustering
parameter$nrCluster <- 10
parameter$pthreshold <- 0.2


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

mdataf <- cp_filterforNA(mdata)
RESULTS$dataDims <- c(RESULTS$dataDims, nrPortNoNas = nrow(mdataf))

# scale matrix rows
scaledM <- t(scale(t(mdataf),center = parameter$row_center, scale = parameter$row_scale))
# ' some alternatives




if (parameter$clustering == "hclust") {
  resClust <- cp_clusterHClustEuclideanDist(scaledM,nrCluster = parameter$nrCluster)
} else if (parameter$clustering == "hclustdeepsplit") {
  resClust <- cp_clusterHClustEuclideanDistDeepslit(scaledM)
} else if (parameter$clustering == "DPA") {
  resClust <- cp_clusterDPAEuclideanDist(scaledM)
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
outfile <- paste0(parameter$projectID,"_",
                  as.character(parameter$workunitID),"_",
                  as.character(outfile),"_",
                  parameter$clustering)


saveRDS(RESULTS,
        file = file.path(parameter$outpath, paste0(outfile,".Rds")))


## create summary.

# input zip organizm outpath

# zipfilename_clusteringname.txt
## 1 row per zip
# Zipfile,
# name of the clustering algorithm
# nr. proteins,
# nr. of samples,
# nr. of protein clusters,
# NR of GS with padj < 0.25
# NR of GS with padj < 0.1
# median CV coefficient of variation for raw data
# median SD after normalizations

# TODO in output 2 and 3.
#ClusterId use integers.

output2 <- lapply(RESULTS$resGOEnrich,
                  function(x){res <- as.data.frame(x$clustProf); res$GS <- x$mt; res})
output2 <- dplyr::bind_rows(output2)

output2 <- tibble::add_column(output2,
                              projectID = parameter$projectID,
                              workunitID = parameter$workunitID,
                              zipfile = basename(parameter$inputMQfile),
                              clustering  = parameter$clustering, .before = 1 )

output2 <- output2 %>% filter(p.adjust < parameter$pthreshold)


RESULTS$nr.of.GS.025 <- output2 %>% filter(p.adjust < 0.25) %>% nrow
RESULTS$nr.of.GS.01 <- output2 %>% filter(p.adjust < 0.1) %>% nrow

# output :
output1 <- data.frame(
  projectID = parameter$projectID,
  workunitID = parameter$workunitID,
  zipfile = basename(parameter$inputMQfile),
  clustering  = parameter$clustering,
  projectID = parameter$projectID,
  workunitID = parameter$workunitID,
  nr.proteins = RESULTS$dataDims["nrPort"],
  nr.proteins.NA.filtered = RESULTS$dataDims["nrPortNoNas"],
  nr.UniprotIDs = RESULTS$dataDims["UniprotExtract"],
  nr.ENTREZIDS = RESULTS$dataDims["ENTREZGENEID"],
  nr.samples = nrow(RESULTS$prot$factors()),
  nr.of.clusters = RESULTS$nrCluster,
  nr.of.GS.025 = RESULTS$nr.of.GS.025,
  nr.of.GS.01 = RESULTS$nr.of.GS.01,
  median.CV = RESULTS$CV.50,
  median.sd = RESULTS$SD.50,
  id.mapping.service = RESULTS$id.mapping.service
)

readr::write_tsv(output1, file = file.path(parameter$outpath, paste0("Summary_", outfile, '.tsv')))
# GS_zipfilename_clusteringname.txt

readr::write_tsv(output2 , file = file.path(parameter$outpath, paste0("GS_", outfile, ".tsv")))

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

tmp <- RESULTS$prot$to_wide()$data
output3 <- right_join(output3, tmp, by = "protein_Id")

protfilename <- paste0("Protein_" , outfile, '.tsv')
readr::write_tsv(output3, file = file.path(parameter$outpath, protfilename))

rmarkdown::render("profileClusters_V2.Rmd",
                  params = list(resultsxx = RESULTS, parametersxx = parameter))

file.copy("profileClusters_V2.html",
          file.path(parameter$outpath, paste0("HTML_",outfile, ".html")),overwrite = TRUE)


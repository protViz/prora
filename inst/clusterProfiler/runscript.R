library(tidyverse)
library(prolfqua)
#library(dendextend)
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
} else if ( length(args) == 1) {
  # read yaml and extract
} else {
  print("::::USING COMMAND LINE ARGS:::")
  parameter$inputMQfile <- args[1]
  parameter$organism <- args[2]
  parameter$outpath <- args[3]
  parameter$clustering <- args[4]
  print(parameter)
}



parameter$clustering <- "hclust"

# related to data preprocessing
parameter$peptide  = FALSE
parameter$transform = "log2scaled" # "log2"

parameter$row_center <- TRUE
parameter$row_scale <- TRUE

# related to clustering
parameter$nrCluster <- 10
parameter$pthreshold <- 0.2


# results will be saved as rds file
results <- list()
if (!dir.exists(parameter$outpath)) {
  dir.create(parameter$outpath)
}

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
  tmp <- prolfqua::tidyMQ_Peptides(parameter$inputMQfile)
  config <- prolfqua::create_config_MQ_peptide()
} else {
  tmp <- prolfqua::tidyMQ_ProteinGroups(parameter$inputMQfile)
  tmp <- tmp %>% filter(.data$nr.peptides > 1)
  atable <- AnalysisTableAnnotation$new()
  atable$fileName = "raw.file"
  atable$hierarchy[["protein_Id"]] <- c("majority.protein.ids")
  atable$hierarchyDepth <- 1
  atable$setWorkIntensity("mq.protein.intensity")
  anaparam <- AnalysisParameters$new()
  config <- AnalysisConfiguration$new(atable, anaparam)
}

tmp$AllSamples <- "AllSamples"
tmp$file <- tmp$raw.file
config$table$factors[["all"]] = "AllSamples"
config$table$factors[["file"]] = "file"

sdata <- setup_analysis(tmp, config)
lfqd <- LFQData$new(sdata, config)
lfqd$filter_proteins_by_peptide_count()

lfqd$remove_small_intensities()
lfqd$hierarchy_counts()

# different ways to normalize data
if (parameter$peptide & parameter$transform == "log2") {
  print("ONLY log2 transform intensities.")
  tr <- lfqd$get_Transformer()
  tr$intensity_array(log2)
  transformed <- tr$lfq
  ag <- transformed$get_Aggregator()
  ag$medpolish()
  prot <- ag$lfq_agg

  tmp <- prot$get_Stats()$stats()
  results$SD.50 <- median(tmp$sd, na.rm = TRUE)

} else if (parameter$peptide & parameter$transform == "log2scaled") {
  tr <- lfqd$get_Transformer()
  transformed <- tr$log2_robscale()
  ag <- transformed$get_Aggregator()
  ag$medpolish()
  prot <- ag$lfq_agg
  tr <- prot$get_Transformer()
  tr$intensity_matrix(.func = prolfqua::robust_scale)
  prot <- tr$lfq

  tmp <- prot$get_Stats()$stats()
  results$SD.50 <- median(tmp$sd, na.rm = TRUE)

} else if (!parameter$peptide & parameter$transform == "log2") {
  st <- lfqd$get_Stats()
  tmp <- st$stats()
  results$CV.50 <- median(tmp$CV, na.rm = TRUE)


  tr <- lfqd$get_Transformer()
  tr$intensity_array(log2)
  prot <- tr$lfq
  tmp <- prot$get_Stats()$stats()

  results$SD.50 <- median(tmp$sd, na.rm = TRUE)

} else if (!parameter$peptide & parameter$transform == "log2scaled") {
  st <- lfqd$get_Stats()
  tmp <- st$stats()
  results$CV.50 <- median(tmp$CV, na.rm = TRUE)

  tr <- lfqd$get_Transformer()
  prot <- tr$log2_robscale()
  tmp <- prot$get_Stats()$stats()

  results$SD.50 <- median(tmp$sd, na.rm = TRUE)
}


results$prot <- prot

wide <- prot$to_wide(as.matrix = TRUE)
mdata <- wide$data
dim(mdata)

### Cluster ----
# what are possible input/outputs?
# dendrogram
clusterHClustEuclideanDist <- function(x, nrCluster , method = "complete"){
  distJK <- prora::dist_JK(x)
  bb <- hclust(distJK,method = method)
  k <- cutree(bb, k = nrCluster)
  dend <- as.dendrogram(bb)
  clusterAssignment <- data.frame(protein_Id = rownames(x), Cluster =  k)
  return(list(dendrogram = dend, clusterAssignment = clusterAssignment))
}

# remove porteins with more NA's than in 60% of samples
filterforNA <- function(mdata){
  na <- apply(mdata, 1, function(x){sum(is.na(x))})
  nc <- ncol(mdata)
  mdata <- mdata[na < floor(0.5 * nc),]
  return(mdata)
}

results$dataDims <- c(nrPort = nrow(mdata))
mdataf <- filterforNA(mdata)
results$dataDims <- c(results$dataDims, nrPortNoNas = nrow(mdataf))
# scale matrix rows
scaledM <- t(scale(t(mdataf),center = parameter$row_center, scale = parameter$row_scale))
# ' some alternatives



if (parameter$clustering == "hclust") {
  resClust <- clusterHClustEuclideanDist(scaledM,nrCluster = parameter$nrCluster)
}else{
  resClust <- clusterHClustEuclideanDist(scaledM,nrCluster = parameter$nrCluster)
}


results$scaledM <- scaledM
results$dendrogram <- resClust$dendrogram



clusterAssignment <- prora::get_UniprotID_from_fasta_header(resClust$clusterAssignment, idcolumn = "protein_Id")
results$dataDims <- c(results$dataDims,  UniprotExtract = sum(!is.na(clusterAssignment$UniprotID)))

# map to entriz id's
.ehandler = function(e){
  warning("WARN :", e)
  # return string here
  as.character(e)
}

#undebug(prora::map_ids_annotationHub)

clusterAssignment <- tryCatch(prora::map_ids_uniprot(clusterAssignment), error = .ehandler)
results$id.mapping.service <- "UniProt"
if (is.character(clusterAssignment)) {
  clusterAssignment <- tryCatch(prora::map_ids_annotationHub(clusterAssignment, species = parameter$species), error = .ehandler)
  results$id.mapping.service <- "AnnotationHub"
}



results$dataDims <- c(results$dataDims,  ENTREZGENEID = sum(!is.na(clusterAssignment$P_ENTREZGENEID)))


# clusterProfiler ----
#clusterAssignmentF
results$clusterAssignment <- clusterAssignment
clusterB <- clusterAssignment %>% filter(!is.na(P_ENTREZGENEID))
# check mapping efficiency.

#clusterB$clusterID <- paste0("Cluster", clusterB$clusterID)
clusterProfilerinput <- split(clusterB$P_ENTREZGENEID, clusterB$Cluster)
length(clusterProfilerinput)


compClust <- function(clusterProfilerinput,
                      orgDB,
                      universe, ont="BP" ){
  clustProf = tryCatch(clusterProfiler::compareCluster(
    clusterProfilerinput,
    fun = "enrichGO",
    OrgDb =  orgDB,
    ont = ont,
    universe = universe,
    pvalueCutoff = 1,
    qvalueCutoff = 0.2,
    pAdjustMethod = "BH"), error = .ehandler)
  return(clustProf)

}

resGOEnrich <- list()
mt <- "GO Biological Process"
clustProf <- compClust(clusterProfilerinput,orgDB,clusterB$P_ENTREZGENEID,ont = "BP")
if (!is.character(clustProf)) {
  resGOEnrich$BP <- list(mt = mt,clustProf = clustProf)
}

mt <- "GO Molecular Function"
clustProf <- compClust(clusterProfilerinput, orgDB, clusterB$P_ENTREZGENEID,ont = "MF")
if (!is.character(clustProf)) {
  resGOEnrich$MF <- list(mt = mt, clustProf = clustProf)
}

mt <- "GO Cellular Component"
clustProf <- compClust(clusterProfilerinput, orgDB, clusterB$P_ENTREZGENEID,ont = "CC")
if (!is.character(clustProf)) {
resGOEnrich$CC <- list(mt = mt, clustProf = clustProf)
}

results$resGOEnrich <- resGOEnrich


outfile <- tools::file_path_sans_ext(basename(parameter$inputMQfile))
saveRDS(results,
        file = file.path(parameter$outpath, paste0(outfile, "_", parameter$clustering,".Rds")))



rmarkdown::render("profileClusters_V2.Rmd",
                  params = list(resultsxx = results))
file.copy("profileClusters_V2.html",
          file.path(parameter$outpath, paste0(outfile,"_", parameter$clustering, ".html")),overwrite = TRUE)

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

output2 <- lapply(results$resGOEnrich,
                  function(x){res <- as.data.frame(x$clustProf); res$GS <- x$mt; res})
output2 <- bind_rows(output2)

output2 <- tibble::add_column(output2,zipfile = basename(parameter$inputMQfile),
                              clustering  = parameter$clustering, .before = 1 )

output2 <- output2 %>% filter(p.adjust < parameter$pthreshold)
colnames(output2)
nr.of.GS.025 <- output2 %>% filter(p.adjust < 0.25) %>% nrow
nr.of.GS.01 <- output2 %>% filter(p.adjust < 0.1) %>% nrow



# output :
output1 <- data.frame(
  zipfile = basename(parameter$inputMQfile),
  clustering  = parameter$clustering,
  nr.proteins = results$dataDims["nrPort"],
  nr.proteins.NA.filtered = results$dataDims["nrPortNoNas"],
  nr.UniprotIDs = results$dataDims["UniprotExtract"],
  nr.ENTREZIDS = results$dataDims["ENTREZGENEID"],
  nr.samples = nrow(results$prot$factors()),
  nr.of.clusters = parameter$nrCluster,
  nr.of.GS.025 = nr.of.GS.025,
  nr.of.GS.01 = nr.of.GS.01,
  median.CV = results$CV.50,
  median.sd = results$SD.50,
  id.mapping.service = results$id.mapping.service
)

write_tsv(output1, file = file.path(parameter$outpath, paste0("Summary", outfile, "_", parameter$clustering, '.tsv')))

# GS_zipfilename_clusteringname.txt

gsfilename <- paste0("GS_" , outfile, "_", parameter$clustering, '.tsv')
readr::write_tsv(output2 , file = file.path(parameter$outpath, gsfilename))

# Protein_zipfilename_clusteringname.txt
# 300 - 5000 rows
# zipfile
# clusterinng
# Protein IDS
# cluster assignments
# Protein intensities sample_Id

protfilename <- paste0("Protein_" , outfile, "_", parameter$clustering, '.tsv')
output3 <- results$clusterAssignment
tmp <- results$prot$to_wide()$data
output3 <- right_join(output3, tmp, by = "protein_Id")
readr::write_tsv(output3, file = file.path(parameter$outpath, protfilename))


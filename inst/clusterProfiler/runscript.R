library(tidyverse)
library(prolfqua)
library(dendextend)
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(prora)

datadir <- file.path(find.package("prolfquaData") , "quantdata")



# parameters
parameter <- list()
parameter$inputMQfile <-  file.path(datadir, "MAXQuant_ComboCourse_p2691_March_2018_WU183012.zip")
parameter$organism <- "yeast" # "human", "mouse"
parameter$outpath = "."


# related to data preprocessing
parameter$peptide  = FALSE
parameter$transform = "log2scaled" # "log2"

parameter$row_center <- TRUE
parameter$row_scale <- TRUE

# related to clustering
parameter$nrCluster <- 10



# resutls will be saved as rds file
results <- list()
dir.create(parameter$outpath)

if (parameter$organism == "yeast") {
  orgDB <- "org.Sc.sgd.db"
} else if (parameter$organism == "human") {
  orgDB <- "org.Hs.eg.db"
} else if (parameter$organism == "mouse") {
  orgDB <- "org.Mm.eg.db"
}


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


tmp$file <- tmp$raw.file
config$table$factors[["file"]] = "file"

sdata <- setup_analysis(tmp, config)
lfqd <- LFQData$new(sdata, config)
lfqd$filter_proteins_by_peptide_count()

lfqd$hierarchy_counts()
lfqd$remove_small_intensities()

# different ways to normalize data
if (parameter$peptide & parameter$transform == "log2") {
  print("ONLY log2 transform intensities.")
  tr <- lfqd$get_Transformer()
  tr$intensity_array(log2)
  transformed <- tr$lfq
  ag <- transformed$get_Aggregator()
  ag$medpolish()
  prot <- ag$lfq_agg
} else if (parameter$peptide & parameter$transform == "log2scaled") {
  tr <- lfqd$get_Transformer()
  transformed <- tr$log2_robscale()
  ag <- transformed$get_Aggregator()
  ag$medpolish()
  prot <- ag$lfq_agg
  tr <- prot$get_Transformer()
  tr$intensity_matrix(.func = prolfqua::robust_scale)
  prot <- tr$lfq
} else if (!parameter$peptide & parameter$transform == "log2") {
  tr <- lfqd$get_Transformer()
  tr$intensity_array(log2)
  prot <- tr$lfq

} else if (!parameter$peptide & parameter$transform == "log2scaled") {
  tr <- lfqd$get_Transformer()
  prot <- tr$log2_robscale()
}


results$prot <- prot

wide <- prot$to_wide(as.matrix = TRUE)
mdata <- wide$data
dim(mdata)

### Cluster ----
# what are possible input/outputs?
# dendrogram
clusterHClustEuclideanDist <- function(mdata, nrCluster ){

  distJK <- prora::dist_JK(scaledM)
  bb <- hclust(distJK,method = )
  k <- cutree(bb, k = nrCluster)
  dend <- as.dendrogram(bb)
  k2 <- cutree(dend, k = parameter$nrCluster)
  all.equal(k, k2)

  clusterAssignment <- data.frame(protein_Id = rownames(scaledM), clusterID =  k)
  return(list(dendrogram = dend, clusterAssignment = clusterAssignment))
}

na <- apply(mdata, 1, function(x){sum(is.na(x))})
nc <- ncol(mdata)
mdata <- mdata[na < floor(0.6 * nc),]
scaledM <- t(scale(t(mdata),center = parameter$row_center, scale = parameter$row_scale))

resClust <- clusterHClustEuclideanDist(scaledM,nrCluster = parameter$nrCluster)

clusterAssignment <- prora::get_UniprotID_from_fasta_header(resClust$clusterAssignment, idcolumn = "protein_Id")
clusterAssignment <- prora::map_ids_uniprot(clusterAssignment)

clusterAssignment %>% filter(is.na(P_ENTREZGENEID))

results$clusterAssignment <- clusterAssignment
results$dendrogram <- resClust$dendrogram

# clusterProfiler ----
clusterB <- na.omit(clusterAssignment)
# check mapping efficiency.

clusterB$clusterID <- paste0("Cluster", clusterB$clusterID)
clusterProfilerinput <- split(clusterB$P_ENTREZGENEID, clusterB$clusterID)

resGOEnrich <- list()
mt <- "GO Biological Process"
resGOEnrich$BP <- list(mt = mt,
                       clustProf = clusterProfiler::compareCluster(
                         clusterProfilerinput,
                         fun = "enrichGO",
                         OrgDb =  orgDB,
                         ont = "BP",
                         universe = clusterB$P_ENTREZGENEID,
                         pvalueCutoff = 1,
                         qvalueCutoff = 0.2,
                         pAdjustMethod = "BH"))

mt <- "GO Molecular Function"
resGOEnrich$MF <- list(mt = mt, clustProf = clusterProfiler::compareCluster(
  clusterProfilerinput,
  fun = "enrichGO",
  OrgDb =  orgDB,
  ont = "MF",
  universe = clusterB$P_ENTREZGENEID,
  pvalueCutoff = 1,
  qvalueCutoff = 0.2,
  pAdjustMethod = "BH"))


mt <- "GO Cellular Component"
resGOEnrich$CC <- list(mt = mt, clustProf = clusterProfiler::compareCluster(
  clusterProfilerinput,
  fun = "enrichGO",
  OrgDb =  orgDB,
  ont = "CC",
  universe = clusterB$P_ENTREZGENEID,
  pvalueCutoff = 1,
  qvalueCutoff = 0.2,
  pAdjustMethod = "BH"))

results$resGOEnrich <- resGOEnrich


outfile <- tools::file_path_sans_ext(basename(parameter$inputMQfile))
saveRDS(results, file = file.path(parameter$outpath, paste0(outfile,".Rds")))


rmarkdown::render("profileClusters_V2.Rmd", params = list(resultsxx = results))


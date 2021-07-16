### Cluster ----
# what are possible input/outputs?
# dendrogram

#' cluster using hclust
#' @export
#'
cp_clusterHClustEuclideanDist <- function(x, nrCluster , method = "complete"){
  distJK <- prora::dist_JK(x)
  bb <- hclust(distJK,method = method)
  k <- cutree(bb, k = nrCluster)
  dend <- as.dendrogram(bb)
  clusterAssignment <- data.frame(protein_Id = rownames(x), Cluster =  k)
  return(list(dendrogram = dend, clusterAssignment = clusterAssignment, nrCluster = nrCluster))
}

#' cluster using hclust and deepsplit
#' @export
#'
cp_clusterHClustEuclideanDistDeepslit <- function(x,  method = "complete"){
  distJK <- prora::dist_JK(x)
  bb <- hclust(distJK,method = method)
  #cat("MIN cluster size ", min(50, nrow(x)/10))
  k <- dynamicTreeCut::cutreeDynamic(
    bb,
    method = "hybrid",
    deepSplit = FALSE,
    distM = as.matrix(distJK),
    minClusterSize = min(50, nrow(x)/10) , verbose = 0)

  dend <- as.dendrogram(bb)
  clusterAssignment <- data.frame(protein_Id = rownames(x), Cluster =  k)
  return(list(dendrogram = dend, clusterAssignment = clusterAssignment, nrCluster = max(k)))
}

#' Cluster using DPA
#' @export
#'
cp_clusterDPAEuclideanDist <- function(mdata, Z = 1){
  distJK <- prora::dist_JK(mdata)
  DPAresult <- runDPAclustering(as.matrix(distJK), Z = Z)
  k <- DPAresult$labels
  maxD <- max(DPAresult$density)
  topography <- DPAresult$topography
  bb <- DPAclustR::plot_dendrogram(k,
                                   topography,
                                   maxD,
                                   popmin = 0,
                                   method = "average")
  dend <- as.dendrogram(bb)
  clusterAssignment <- data.frame(protein_Id = rownames(mdata), Cluster =  k)
  return(list(dendrogram = dend, clusterAssignment = clusterAssignment, nrCluster = max(k)))
}

#' remove proteins with more NA's than in 60% of samples
#' @export
#'
cp_filterforNA <- function(mdata){
  na <- apply(mdata, 1, function(x){sum(is.na(x))})
  nc <- ncol(mdata)
  mdata <- mdata[na < floor(0.5 * nc),]
  return(mdata)
}

#' compute clusters
#' @export
#'
cp_compClust <- function(clusterProfilerinput,
                         orgDB,
                         universe, ont="BP" ){
  # map to entriz id's
  .ehandler = function(e) {
    warning("WARN :", e)
    # return string here
    as.character(e)
  }

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


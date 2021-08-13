### Cluster ----
# what are possible input/outputs?
# dendrogram

#' cluster using hclust
#' @export
#'
cp_clusterHClustEuclideanDist <- function(x, nrCluster , method = "complete", JK = TRUE){
  distJK <- if (JK) {
     prora::dist_JK(x)
  } else {
    dist(x)
  }

  bb <- hclust(distJK, method = method)
  k <- cutree(bb, k = nrCluster)
  dend <- as.dendrogram(bb)
  clusterAssignment <- data.frame(protein_Id = rownames(x), Cluster =  k)
  return(list(dendrogram = dend, clusterAssignment = clusterAssignment, nrCluster = nrCluster))
}

#' cluster using hclust and deepsplit
#' @export
#'

cp_clusterHClustEuclideanDistDeepslit <- function(x,  method = "complete", JK = TRUE){
  if (nrow(x) <= 2) {
    dend <- NULL
    nrCluster <- 1
    clusterAssignment <- data.frame( protein_Id  = as.character(rownames(x)), Cluster = rep(1, nrow(x) ))
  } else {
    distJK <- if (JK) {
      prora::dist_JK(x)
    } else {
      dist(x)
    }

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
    nrCluster <-  max(k)
  }
  return(list(dendrogram = dend, clusterAssignment = clusterAssignment, nrCluster = nrCluster))
}

#' Cluster using DPA
#'
#' @param mdata, input data
#' @param Z, the number of standard deviations fixing the level of
#'            statistical confidence at which one decides to consider
#'            a cluster meaningful. Default value is set to 1.
#' @param metric, string, valid metrics are "precomputed", "euclidean".
#'            If metric is "precomputed", data is assumed to be a distance matrix.
#'            Default is "euclidean".
#' @export
#'
cp_clusterDPAEuclideanDist <- function(mdata, Z = 1, metric = "euclidean", JK = TRUE){
  distJK <- if (JK) {
    prora::dist_JK( mdata )
  } else {
    dist(mdata)
  }
  DPAresult <- DPAclustR::runDPAclustering(as.matrix(distJK), Z = Z, metric = metric)
  k <- DPAresult$labels
  maxD <- max(DPAresult$density)
  topography <- DPAresult$topography
  if ( nrow(topography) > 0 ) {
    bb <- DPAclustR::plot_dendrogram(k,
                                     topography,
                                     maxD,
                                     popmin = 0,
                                     method = "average")

    dend <- as.dendrogram(bb)
  } else {
    dend <- NULL
  }
  clusterAssignment <- data.frame(protein_Id = rownames(mdata), Cluster =  k)
  return( list(dendrogram = dend, clusterAssignment = clusterAssignment, nrCluster = max(k)) )
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


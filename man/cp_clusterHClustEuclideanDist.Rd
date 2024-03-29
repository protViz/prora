% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cphelper.R
\name{cp_clusterHClustEuclideanDist}
\alias{cp_clusterHClustEuclideanDist}
\title{cluster using hclust}
\usage{
cp_clusterHClustEuclideanDist(x, nrCluster, method = "complete", JK = TRUE)
}
\arguments{
\item{x}{data matrix}

\item{nrCluster}{number of clusters}

\item{method}{default complete}

\item{JK}{jack knife resmapling default TRUE}
}
\description{
cluster using hclust
}
\seealso{
Other clusterProfiler: 
\code{\link{cp_addurls_GO}()},
\code{\link{cp_clusterDPAEuclideanDist}()},
\code{\link{cp_clusterHClustEuclideanDistDeepslit}()},
\code{\link{cp_compClust}()},
\code{\link{cp_filterforNA}()}
}
\concept{clusterProfiler}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fgsea_utils.R
\name{fgsea_rank_contrasts}
\alias{fgsea_rank_contrasts}
\title{convert data frame into list of ranks}
\usage{
fgsea_rank_contrasts(df, ids = "X1", score = "X2", contrast = "contrast")
}
\arguments{
\item{df}{data frame}

\item{ids}{column with ids}

\item{score}{column with scores}

\item{contrast}{column with contrast name.}
}
\description{
convert data frame into list of ranks
}
\seealso{
Other fgsea: 
\code{\link{fgsea_leading_edge_too_char}()},
\code{\link{fgsea_msigdb_collections}()},
\code{\link{fgsea_msigdb}()},
\code{\link{fgsea_rank}()},
\code{\link{getMsigdbGenesets}()},
\code{\link{run_fgsea_for_allContrasts}()},
\code{\link{run_fgsea_for_allGeneSets}()}
}
\concept{fgsea}

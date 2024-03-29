% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigoraWrappR.R
\name{sigoraWrappR}
\alias{sigoraWrappR}
\title{sigora Wrapper function}
\usage{
sigoraWrappR(
  data,
  threshold = 0.5,
  score_col = "",
  GPSrepos = sigora::load_data("kegH"),
  db = "",
  greater_than = TRUE,
  idmap = sigora::idmap
)
}
\arguments{
\item{data}{input data.frame (at least two columns, first column containing IDs,
other columns numerical ranks, i.e. fold changes)}

\item{threshold}{fold change threshold above which (in absolute terms) a protein
is considered differentially regulated}

\item{score_col}{Name of the fold change column, in case the input file contains
multiple contrasts}

\item{GPSrepos}{GPS repository used as background, can be
generated via \code{\link{sigoraWrappR}}}

\item{db}{database used for the generation of the GPS repository}

\item{greater_than}{Logical. Whether fc_threshold should be applied as
greater_than (default is \code{TRUE}) or strictly less than (\code{FALSE})}

\item{idmap}{id mapping table, dafault sigora::idmap}
}
\value{
Returns a \code{list} containing the following elements:
\itemize{
   \item \code{sigora}: \code{data.frame} containing sigora results
   \item \code{ora}: \code{data.frame} containing ora results
   \item \code{fc_threshold}: fc_threshold used
   \item \code{GPS_repository}: GPS repository used, to be reported
   \item \code{database}: Database used for generating the GPS repository
   \item \code{data}: \code{data.frame} containing Uniprot IDs and fold changes
   \item \code{proteinsAfterFiltering}: numeric, number of proteins after fc filtering
}
}
\description{
Provides a wrapper for \code{sigora::sigora()} combined with \code{sigora::ora()} and
returns a \code{list} used for the generation of an \code{.Rmd} report.
}
\examples{

library(prora)
df <- get_UniprotID_from_fasta_header(prora::exampleContrastData)
myGPSrepo <- makeGPS_wrappR(ids = df$UniprotID,target = "KEGG")
res <- sigoraWrappR(df,score_col = "estimate", GPSrepos = myGPSrepo$gps,
                     threshold = 0.5)

myGPSrepoGO <- makeGPS_wrappR(df$UniprotID, target = "GO")
res <- sigoraWrappR(df,score_col = "estimate", GPSrepos = myGPSrepoGO$gps,
                     threshold = 0.5)
\dontrun{
  res$ora
res$sigora$summary_results
res$data
GPSrepo <-  myGPSrepoGO$gps
}
}

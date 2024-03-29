% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runWebGestalt.R
\name{runWebGestaltORA}
\alias{runWebGestaltORA}
\title{workflow for ORA}
\usage{
runWebGestaltORA(
  data,
  fpath,
  organism = "hsapiens",
  ID_col = "UniprotID",
  target = "geneontology_Biological_Process",
  threshold = 0.59,
  greater = TRUE,
  nperm = 10,
  score_col = "estimate",
  outdir = "WebGestalt_ORA",
  subdir_name = NULL,
  interestGeneType = "uniprotswissprot",
  contrast_name = fpath
)
}
\arguments{
\item{data}{data}

\item{fpath}{file path to write to}

\item{organism}{organism}

\item{ID_col}{column name containing IDs}

\item{target}{target database, default: "geneontology_Biological_Process"}

\item{threshold}{threshold for estimate}

\item{greater}{indicating direction of threshold}

\item{nperm}{number of permutations}

\item{score_col}{column name containing estimates}

\item{outdir}{output directory}

\item{subdir_name}{subdirectory name}

\item{interestGeneType}{what type of identifier default : "uniprotswissprot"}

\item{contrast_name}{default fpath}
}
\description{
workflow for ORA
}

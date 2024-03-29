% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msigDB_helper.R
\name{getMsigdbGenesets}
\alias{getMsigdbGenesets}
\title{getMsigdbGenesets}
\usage{
getMsigdbGenesets(msigCollection, species)
}
\arguments{
\item{msigCollection}{data.frema with columns gs_cat, gs_subcat}

\item{species}{e.g. "Homo sapiens"}
}
\description{
getMsigdbGenesets
}
\examples{
library(tidyverse)
library(msigdbr)
msigdbr::msigdbr_species()
species <- "Homo sapiens"
species <- "Mus musculus"

hallmark <- msigdbr_collections() \%>\% filter(.data$gs_cat == "H")

#hallmark$gs_subcat <- "HALLMARK"
C5 <- bind_rows( {msigdbr_collections() \%>\%
 filter(.data$gs_cat == "C5") \%>\%
 filter(grepl("^GO:", .data$gs_subcat))},
 hallmark,
 {msigdbr_collections() \%>\% filter(.data$gs_subcat == "CP:KEGG")} )

C5
fgseaGSlist <- prora::getMsigdbGenesets(C5, species)

}
\seealso{
Other fgsea: 
\code{\link{fgsea_leading_edge_too_char}()},
\code{\link{fgsea_msigdb_collections}()},
\code{\link{fgsea_msigdb}()},
\code{\link{fgsea_rank_contrasts}()},
\code{\link{fgsea_rank}()},
\code{\link{run_fgsea_for_allContrasts}()},
\code{\link{run_fgsea_for_allGeneSets}()}
}
\concept{fgsea}

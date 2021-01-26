# Test sigora

library(WebGestaltR)
library(tidyverse)
library(org.Hs.eg.db)
library(sigora)
library(GO.db)
library(slam)
library(prora)

dd <- prora::exampleContrastData
dfiles <- "example_data.txt"

# Parameters --------------------------------------------------------------

organism <- "hsapiens"
ID_col <- "UniprotID"
fc_col <- "estimate"
target_SIGORA <- c("GO", "KEGG", "reactome")
map_col <- "GO"
fc_threshold <- 0.5
greater <- TRUE


# Run ---------------------------------------------------------------------

fpath_se <- tools::file_path_sans_ext(basename(dfiles))
odir <- make.names(fpath_se)

if (!dir.exists(odir)) {
  dir.create(odir)
}

colnames(dd) <- make.names(colnames(dd))
dd <- get_UniprotID_from_fasta_header(dd) %>%
  filter(!is.na(UniprotID))
sapply(target_SIGORA, function(x) {
  prora::runSIGORA(
    data = dd,
    target = x,
    fc_col = fc_col,
    fc_threshold = fc_threshold,
    greater = greater,
    outdir = file.path(odir, "sigORA")
  )
})

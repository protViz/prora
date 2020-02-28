# Test GSEA

library(WebGestaltR)
library(tidyverse)
library(org.Hs.eg.db)
library(sigora)
library(GO.db)
library(slam)
library(fgcz.gsea.ora)

dd <- fgcz.gsea.ora::exampleContrastData
dfiles <- "example_data.txt"

# Parameters --------------------------------------------------------------

organism <- "hsapiens"
ID_col <- "UniprotID"
fc_col <- "estimate"
target_GSEA <- c(
  "geneontology_Biological_Process",
  "geneontology_Cellular_Component",
  "geneontology_Molecular_Function"
)
fc_threshold <- 0.5
greater <- TRUE
nperm <- 10


# Run ---------------------------------------------------------------------

fpath_se <- tools::file_path_sans_ext(basename(dfiles))
odir <- make.names(fpath_se)

if (!dir.exists(odir)) {
  dir.create(odir)
}

colnames(dd) <- make.names(colnames(dd))
dd <- get_UniprotID_from_fasta_header(dd) %>%
  filter(!is.na(UniprotID))

sapply(target_GSEA, function(x) {
  fgcz.gsea.ora::runWebGestaltGSEA(
    data = dd,
    fpath = "",
    ID_col = ID_col,
    fc_col = fc_col,
    organism = organism,
    target = x,
    nperm = nperm,
    outdir = file.path(odir, "GSEA")
  )
})

library(WebGestaltR)
library(tidyverse)
library(org.Hs.eg.db)
library(sigora)
library(GO.db)
library(slam)
library(fgczgseaora)


# Files -------------------------------------------------------------------

#data("exampleContrastData", package = "fgczgseaora")
dd <- fgczgseaora::exampleContrastData
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

target_GSEA <- target_GSEA[1]
map_col <- "GO"
fc_threshold <- 0.5
greater <- TRUE
nperm <- 10


# Run ---------------------------------------------------------------------

fpath_se <- tools::file_path_sans_ext(basename(dfiles))
odir <- paste0(make.names(fpath_se),"_ORA")

if (!dir.exists(odir)) {
  dir.create(odir)
}

colnames(dd) <- make.names(colnames(dd))
filtered_dd <- getUniprotFromFastaHeader(dd) %>%
  filter(!is.na(UniprotID))


sapply(target_GSEA, function(x) {
  fgczgseaora:::.runWebGestaltORA(
    data = filtered_dd,
    fpath = "",
    organism = organism,
    ID_col = ID_col,
    target = x,
    map_col = map_col,
    threshold = fc_threshold,
    greater = greater,
    nperm = nperm,
    fc_col = fc_col,
    outdir = file.path(odir, "WebGestaltORA")
  )
})

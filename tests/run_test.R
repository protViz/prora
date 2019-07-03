library(WebGestaltR)
library(tidyverse)
library(org.Hs.eg.db)
library(sigora)
library(GO.db)
library(slam)
library(fgczgseaora)


# Files -------------------------------------------------------------------

data("exampleContrastData", package = "fgczgseaora")
dd <- exampleContrastData
dfiles <- "example_data.txt"

# Parameters --------------------------------------------------------------

organism <- "hsapiens"
ID_col <- "UniprotID"
fc_col <- "estimate.Age.class..Old...Young"
target_GSEA <- c(
  "geneontology_Biological_Process",
  "geneontology_Cellular_Component",
  "geneontology_Molecular_Function"
)[1]
target_SIGORA <- c("GO", "KEGG", "reactome")[1]
map_col <- "GO"
fc_threshold <- 0.5
greater <- TRUE
nperm <- 10


# Run ---------------------------------------------------------------------

for (i in seq_len(length(dfiles))) {
  fpath_se <- tools::file_path_sans_ext(basename(dfiles[i]))
  odir <- make.names(fpath_se)

  if (!dir.exists(odir)) {
    dir.create(odir)
  }

  colnames(dd) <- make.names(colnames(dd))
  dd <- getUniprotFromFastaHeader(dd) %>%
    filter(!is.na(UniprotID))

  sapply(target_GSEA, function(x) {
    fgczgseaora:::.runGSEA(
      data = dd,
      fpath = "",
      ID_col = ID_col,
      fc_col = fc_col,
      organism = organism,
      target = x,
      map_col = map_col,
      nperm = nperm,
      outdir = file.path(odir, "GSEA")
    )
  })

  sapply(target_GSEA, function(x) {
    fgczgseaora:::.runWebGestaltORA(
      data = dd,
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

  if (organism == "hsapiens")
    sapply(target_SIGORA, function(x) {
      fgczgseaora:::.runSIGORA(
        data = dd,
        target = x,
        fc_col = fc_col,
        fc_threshold = fc_threshold,
        greater = greater,
        outdir = file.path(odir, "sigORA")
      )
    })
}

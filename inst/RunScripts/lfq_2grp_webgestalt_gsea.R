library(WebGestaltR)
library(tidyverse)
library(org.Hs.eg.db)
library(sigora)
library(GO.db)
library(slam)
library(fgczgseaora)
library(readr)

# Files -------------------------------------------------------------------

grp2report <- "data/2Grp_CF_a_vs_CF_b.txt"
result_dir <- "gsea_ora_results"


target_GSEA <- c(
  "geneontology_Biological_Process",
  "geneontology_Cellular_Component",
  "geneontology_Molecular_Function"
)


#target_SIGORA <- target_SIGORA[1]
# Parameters --------------------------------------------------------------


ID_col <- "TopProteinName"
fc_col <- "log2FC"
nperm <- 100

#
fpath_se <- tools::file_path_sans_ext(basename(grp2report))
odir <- file.path(result_dir , make.names(fpath_se))
if (!dir.exists(odir)) {
  dir.create(odir)
}


dd <- read_tsv(grp2report)
dd <- dd %>% select_at(c(ID_col, fc_col))


filtered_dd <- getUniprotFromFastaHeader(dd,idcolumn = ID_col) %>%
  filter(!is.na(UniprotID))

filtered_dd <- na.omit(filtered_dd)


res <- lapply(target_GSEA, function(x) {
  message(x)
  fgczgseaora:::.runGSEA(
    data = filtered_dd,
    fpath = "",
    ID_col = "UniprotID",
    fc_col = fc_col,
    organism = organism,
    target = x,
    nperm = nperm,
    outdir = file.path(odir, "GSEA")
  )
})



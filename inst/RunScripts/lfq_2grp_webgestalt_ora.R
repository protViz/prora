#lfq_2grp_webgestalt_ORA.R
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


ID_col <- "TopProteinName"
fc_col <- "log2FC"
nperm <- 10

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
  message("\n",x,"\n")
  fgczgseaora:::.runWebGestaltORA(

    data = filtered_dd,
    fpath = "",
    organism = organism,
    ID_col = "UniprotID",
    target = x,
    threshold = fc_threshold,
    greater = greater,
    nperm = nperm,
    fc_col = fc_col,
    outdir = file.path(odir, "WebGestaltORA")
  )
})



#!/usr/bin/env Rscript

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
args = commandArgs(trailingOnly=TRUE)
if (length(args) < 2) {
  stop("At least two argument must be supplied: grp2file and output_directory", call.=FALSE)
}


grp2report <- args[1]
result_dir <- args[2]

nperm <- 10

if (length(args)==3) {
  # default output file
  nperm <- args[3]
}

ID_col <- "TopProteinName"
fc_col <- "log2FC"
if (length(args) == 5) {
  # default output file
  ID_col <- args[4]
  fc_col <- args[5]
}



target_GSEA <- c(
  "geneontology_Biological_Process",
  "geneontology_Cellular_Component",
  "geneontology_Molecular_Function"
)



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



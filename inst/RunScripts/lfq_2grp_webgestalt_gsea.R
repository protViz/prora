#!/usr/bin/Rscript
"WebGestaltR GSEA

Usage:
  test.R <grp2file> [--organism=<organism>] [--outdir=<outdir>]
        [--nperm=<nperm>] [--ID_col=<ID_col>] [--fc_col=<fc_col>]

Options:
  -o --organism=<organism> organism [default: hsapiens]
  -r --outdir=<outdir> output directory [default: results]
  -n --nperm=<nperm> number of permutations to calculate enrichment scores [default: 10]
  -i --ID_col=<ID_col> Column containing the UniprotIDs [default: TopProteinName]
  -f --fc_col=<fc_col> Column containing the estimates [default: log2FC]

Arguments:
  grp2file  input file
" -> doc

library(docopt)
opt <- docopt(doc)

options(warn = -1)
suppressMessages(library(WebGestaltR))
suppressMessages(library(tidyverse))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(sigora))
suppressMessages(library(GO.db))
suppressMessages(library(slam))
suppressMessages(library(fgczgseaora))
suppressMessages(library(readr))


# Check command args ------------------------------------------------------

cat("\nParameters used:\n\t grp2report:", grp2report <- opt$grp2file, "\n\t",
    "result_dir:", result_dir <- opt[["--outdir"]], "\n\t",
    "  organism:", org <- opt[["--organism"]], "\n\t",
    "     nperm:", nperm <- as.numeric(opt[["--nperm"]]), "\n\t",
    "    ID_col:", ID_col <- opt[["--ID_col"]],"\n\t",
    "    fc_col:", fc_col <- opt[["--fc_col"]] , "\n\n\n")


target_GSEA <- c(
  "geneontology_Biological_Process",
  "geneontology_Cellular_Component",
  "geneontology_Molecular_Function"
)


# Parameters --------------------------------------------------------------

fpath_se <- tools::file_path_sans_ext(basename(grp2report))
odir <- file.path(result_dir , make.names(fpath_se))

if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

if (!dir.exists(odir)) {
  dir.create(odir)
}

dd <- read_tsv(grp2report)
dd <- dd %>% select_at(c(ID_col, fc_col))

filtered_dd <- getUniprotFromFastaHeader(dd, idcolumn = ID_col) %>%
  filter(!is.na(UniprotID))
filtered_dd <- na.omit(filtered_dd)

res <- lapply(target_GSEA, function(x) {
  message(x)
  fgczgseaora:::.runGSEA(
    data = filtered_dd,
    fpath = "",
    ID_col = "UniprotID",
    fc_col = fc_col,
    organism = org,
    target = x,
    nperm = nperm,
    outdir = file.path(odir, "GSEA")
  )
})



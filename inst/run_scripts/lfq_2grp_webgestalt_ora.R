#!/usr/bin/Rscript
"WebGestaltR ORA

Usage:
  test.R <grp2file> [--organism=<organism>] [--outdir=<outdir>] [--log2fc=<log2fc>] [--is_greater=<is_greater>] [--nperm=<nperm>] [--ID_col=<ID_col>] [--fc_col=<fc_col>]

Options:
  -o --organism=<organism> organism [default: hsapiens]
  -r --outdir=<outdir> output directory [default: results_ora]
  -t --log2fc=<log2fc> fc threshold [default: 1]
  -g --is_greater=<is_greater> is greater than log2fc [default: TRUE]
  -n --nperm=<nperm> number of permutations to calculate enrichment scores [default: 50]
  -i --ID_col=<ID_col> Column containing the UniprotIDs [default: TopProteinName]
  -f --fc_col=<fc_col> Column containing the estimates [default: log2FC]

Arguments:
  grp2file  input file
" -> doc


suppressMessages(library(WebGestaltR))
suppressMessages(library(tidyverse))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(sigora))
suppressMessages(library(GO.db))
suppressMessages(library(slam))
suppressMessages(library(fgczgseaora))
suppressMessages(library(readr))
library(docopt)


opt <- docopt(doc)

print(opt)

# Check command args ------------------------------------------------------

cat("\nParameters used:\n\t grp2report:", grp2report <- opt$grp2file, "\n\t",
    "result_dir:", result_dir <- opt[["--outdir"]], "\n\t",
    "  organism:", organism <- opt[["--organism"]], "\n\t",
    "    log2fc:", log2fc <- opt[["--log2fc"]], "\n\t",
    "is_greater:", is_greater <- opt[["--is_greater"]], "\n\t",
    "     nperm:", nperm <- as.numeric(opt[["--nperm"]]), "\n\t",
    "    ID_col:", ID_col <- opt[["--ID_col"]],"\n\t",
    "    fc_col:", fc_col <- opt[["--fc_col"]] , "\n\n\n")


target_GSEA <- c(
  "geneontology_Biological_Process",
  "geneontology_Cellular_Component",
  "geneontology_Molecular_Function"
)


organisms <- listOrganism(hostName = "http://www.webgestalt.org/", cache = NULL)

if(! organism %in% organisms){
  stop("organism : ",organism, "is not in the list of available organisms", paste(organisms, collapse=" ") )
}

# Parameters --------------------------------------------------------------

fpath_se <- tools::file_path_sans_ext(basename(grp2report))
#odir <- file.path(result_dir , make.names(fpath_se))

if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

subdir <- file.path(result_dir, paste0("fc_",log2fc,"_is_g_",is_greater))
if(!dir.exists(subdir)){
  dir.create(subdir)
}



fc_estimates <- read_tsv(grp2report)
fc_estimates <- fc_estimates %>% select_at(c(ID_col, fc_col))

print("Selected columns: ")
print(sample_n(fc_estimates, 10))

filtered_dd <- getUniprotFromFastaHeader(fc_estimates, idcolumn = ID_col) %>%
  filter(!is.na(UniprotID))
filtered_dd <- na.omit(filtered_dd)
print("After ID filtering columns: ")
print(sample_n(filtered_dd, 10))




res <- lapply(target_GSEA, function(x) {
  message("\n",x,"\n")
  fgczgseaora::runWebGestaltORA(
    data = filtered_dd,
    fpath = "",
    organism = organism,
    ID_col = "UniprotID",
    target = x,
    threshold = log2fc,
    greater = is_greater,
    nperm = nperm,
    fc_col = fc_col,
    outdir = subdir
  )
})



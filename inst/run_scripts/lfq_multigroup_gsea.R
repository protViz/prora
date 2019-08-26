#!/usr/bin/Rscript
"WebGestaltR GSEA

Usage:
  test.R <grp2file> [--organism=<organism>] [--outdir=<outdir>] [--nperm=<nperm>]

Options:
  -o --organism=<organism> organism [default: hsapiens]
  -r --outdir=<outdir> output directory [default: results_gsea]
  -n --nperm=<nperm> number of permutations to calculate enrichment scores [default: 50]

Arguments:
  grp2file  input file
" -> doc

library(docopt)
opt <- docopt(doc)

options(warn = -1)
suppressMessages( library(WebGestaltR) )
suppressMessages( library(tidyverse) )
suppressMessages( library(org.Hs.eg.db) )
suppressMessages(library(sigora))
suppressMessages(library(GO.db))
suppressMessages(library(slam))
suppressMessages(library(fgczgseaora))
suppressMessages(library(readr))


# Check command args ------------------------------------------------------

cat("\nParameters used:\n\t grp2report:", grp2report <- opt$grp2file, "\n\t",
    "result_dir:", result_dir <- opt[["--outdir"]], "\n\t",
    "  organism:", organism <- opt[["--organism"]], "\n\t",
    "     nperm:", nperm <- as.numeric(opt[["--nperm"]]), "\n\n\n")


target_GSEA <- c(
  "geneontology_Biological_Process",
  "geneontology_Cellular_Component",
  "geneontology_Molecular_Function"
)

ID_col <- "top_protein"
fc_col <- "estimate"
contrast <- "contrast"

organisms <- listOrganism(hostName = "http://www.webgestalt.org/", cache = NULL)

if(! organism %in% organisms){
  stop("organism : " , organism , "is not in the list of available organisms", paste(organisms, collapse=" ") )
}

# Parameters --------------------------------------------------------------


if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

fc_estimates <- readxl::read_xlsx(grp2report)
fc_estimates <- fc_estimates %>% select_at(c(ID_col, fc_col, contrast))

print("Selected columns: ")
print(sample_n(fc_estimates, 10))

filtered_dd <- getUniprotFromFastaHeader(fc_estimates, idcolumn = ID_col) %>%
  filter(!is.na(UniprotID))
filtered_dd <- na.omit(filtered_dd)
print("After ID filtering columns: ")
print(sample_n(filtered_dd, 10))


filtered_dd_list <- base::split(filtered_dd, filtered_dd$contrast)
contr_names <- names(filtered_dd_list)
contr_names <- gsub(" ","", contr_names)
contr_names <- gsub("-","_vs_", contr_names)
contr_names <- make.names(contr_names)

names(filtered_dd_list) <- contr_names

for(name in names(filtered_dd_list)){
  filtered_dd <- filtered_dd_list[[name]]
  cat("\n\n processing contrast :",name, "\n\n")
  res <- lapply(target_GSEA, function(x) {
    message(x)
    fgczgseaora:::.runGSEA(
      data = filtered_dd,
      fpath = name,
      ID_col = "UniprotID",
      fc_col = fc_col,
      organism = organism,
      target = x,
      nperm = nperm,
      outdir = result_dir
    )
  })
}

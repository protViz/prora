#!/usr/bin/Rscript

"WebGestaltR ORA for multigroup reports

Usage:
  test.R <grp2file> [--organism=<organism>] [--outdir=<outdir>] [--nperm=<nperm>] [--log2fc=<log2fc>] [--is_greater=<is_greater>]

Options:
  -o --organism=<organism> organism [default: hsapiens]
  -r --outdir=<outdir> output directory [default: results_ora]
  -n --nperm=<nperm> number of permutations to calculate enrichment scores [default: 50]
  -t --log2fc=<log2fc> fc threshold [default: 1]
  -g --is_greater=<is_greater> is greater than log2fc [default: TRUE]

Arguments:
  grp2file  input file
" -> doc

library(docopt)
opt <- docopt(doc)

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
    "  organism:", organism <- opt[["--organism"]], "\n\t",
    "    log2fc:", log2fc <- as.numeric(opt[["--log2fc"]]), "\n\t",
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

if(! organism %in% organisms) {
  stop("organism : " , organism , "is not in the list of available organisms", paste(organisms, collapse=" ") )
}

# Parameters --------------------------------------------------------------

fpath_se <- tools::file_path_sans_ext(basename(grp2report))
odir <- file.path(result_dir , make.names(fpath_se))

if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

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




log2fc_s <- c(log2fc, -as.numeric(log2fc))

for(log2fc in log2fc_s){
  is_greater <- if(log2fc > 0 ){TRUE}else{FALSE}

  subdir <- file.path(result_dir, paste0("fc_threshold_",abs(log2fc),"_is_greater_",is_greater))
  if(!dir.exists(subdir)){
    message("created directory : ", subdir, "\n\n")
    dir.create(subdir)
  }


  for(name in names(filtered_dd_list)){
    filtered_dd <- filtered_dd_list[[name]]

    cat("\n\n processing contrast :",name, "\n\n")

    res <- lapply(target_GSEA, function(x) {
      message("\n",x,"\n")
      fgczgseaora:::.runWebGestaltORA(
        data = filtered_dd,
        fpath = name,
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
  }
}

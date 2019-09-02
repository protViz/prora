#!/usr/bin/Rscript
"WebGestaltR GSEA for multigroup reports

Usage:
  test.R <grp2file> [--organism=<organism>] [--outdir=<outdir>] [--nperm=<nperm>] [--idtype=<idtype>] [--ID_col=<ID_col>]

Options:
  -o --organism=<organism> organism [default: hsapiens]
  -r --outdir=<outdir> output directory [default: results_gsea]
  -t --idtype=<idtype> type of id used for mapping [default: uniprotswissprot]
  -i --ID_col=<ID_col> Column containing the UniprotIDs [default: top_protein]
  -n --nperm=<nperm> number of permutations to calculate enrichment scores [default: 50]

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
    "  organism:", organism <- opt[["--organism"]], "\n\t",
    "    idtype:", idtype <- opt[["--idtype"]], "\n\t",
    "    ID_col:", idcolumn <- opt[["--ID_col"]], "\n\t",
    "     nperm:", nperm <- as.numeric(opt[["--nperm"]]), "\n\n\n")


target_GSEA <- c(
  "geneontology_Biological_Process",
  "geneontology_Cellular_Component",
  "geneontology_Molecular_Function"
)

ID_col <- idcolumn
fc_col <- "estimate"
contrast <- "contrast"

organisms <- listOrganism(hostName = "http://www.webgestalt.org/", cache = NULL)

if(! organism %in% organisms){
  cat( paste(organisms, collapse="\n") )
  stop("organism : " , organism , "is not in the list of available organisms!" )
}

idtypes <- listIdType(organism = organism, hostName = "http://www.webgestalt.org/", cache = NULL)

if(! idtype %in% idtypes){
  cat(paste(idtypes, collapse="\n"))
  stop("idtype : " , idtypes , "is not in the list of available idtypes!" )
}

# Parameters --------------------------------------------------------------


if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

fc_estimates <- readxl::read_xlsx(grp2report)
fc_estimates <- fc_estimates %>% select_at(c(ID_col, fc_col, contrast))

print("Selected columns: ")
print(sample_n(fc_estimates, 10))


filtered_dd <- na.omit(fc_estimates)
print("After ID filtering columns: ")
print(sample_n(filtered_dd, 10))


filtered_dd_list <- base::split(filtered_dd, filtered_dd$contrast)
contr_names <- names(filtered_dd_list)
contr_names <- gsub(" ","", contr_names)
contr_names <- gsub("-","_vs_", contr_names)
contr_names <- make.names(contr_names)

names(filtered_dd_list) <- contr_names

res <- list()

for(target in target_GSEA)
{
  res_contrast <- list()
  for(name in names(filtered_dd_list))
  {
    filtered_dd <- filtered_dd_list[[name]]
    cat("\n\n processing target : ",target," for contrast : ", name, "\n\n")

    res_contrast[[name]] <- lapply(target_GSEA, function(x) {
      message(x)
      fgczgseaora::runGSEA(
        data = filtered_dd,
        fpath = name,
        ID_col = ID_col,
        fc_col = fc_col,
        organism = organism,
        target = x,
        nperm = nperm,
        outdir = result_dir,
        interestGeneType = idtype
      )
    })
  }
  res[[target]] <- res_contrast
}

saveRDS(res, file="whats_in_GSEA_res.Rda")

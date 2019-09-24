#!/usr/bin/Rscript

options(nwarnings = 10000)

"WebGestaltR GSEA for multigroup reports

Usage:
  test.R <grp2file> [--organism=<organism>] [--outdir=<outdir>] [--idtype=<idtype>] [--ID_col=<ID_col>]  [--nperm=<nperm>] [--score_col=<score_col>] [--contrast=<contrast>]

Options:
  -o --organism=<organism> organism [default: hsapiens]
  -r --outdir=<outdir> output directory [default: results_gsea]
  -t --idtype=<idtype> type of id used for mapping [default: uniprotswissprot]
  -i --ID_col=<ID_col> Column containing the UniprotIDs [default: UniprotID]
  -n --nperm=<nperm> number of permutations to calculate enrichment scores [default: 500]
  -e --score_col=<score_col> column containing fold changes [default: estimate]
  -c --contrast=<contrast> column containing fold changes [default: contrast]

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
    "     nperm:", nperm <- as.numeric(opt[["--nperm"]]), "\n\t",
    "  contrast:", contrast <- opt[["--contrast"]], "\n\t",
    "  score_col:", score_col <- opt[["--score_col"]], "\n\n\n")


target_GSEA <- c(
  "geneontology_Biological_Process",
  "geneontology_Cellular_Component",
  "geneontology_Molecular_Function"
)

ID_col <- idcolumn
fc_col <- score_col

organisms <- listOrganism(hostName = "http://www.webgestalt.org/", cache = NULL)

if(! organism %in% organisms){
  cat("ERROR !\n")
  cat("Organism : " , organism , "is not in the list of available organisms!")
  cat("List of available orginisms\n")
  cat( paste(organisms, collapse="\n") )
  stop("ERROR !\n" )
}

idtypes <- listIdType(organism = organism, hostName = "http://www.webgestalt.org/", cache = NULL)

if(! idtype %in% idtypes){
  cat("ERROR !\n")
  cat("idtype : " , idtype , "is not in the list of available idtypes!\n" )
  cat("list of available idtypes?\n")
  cat(paste(idtypes, collapse="\n"))
  stop("ERROR !\n")
}

# Parameters --------------------------------------------------------------


result_dir <- paste0(result_dir,"_",format(Sys.time(), '%d%b%Y_%H%M%S'))
cat("creating dir ", result_dir,"\n")
if(dir.exists(result_dir)){
  unlink(result_dir, recursive = TRUE)
}
dir.create(result_dir)

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
    cat("\n\n PROCESSING TARGET : ",target," FOR CONTRAST : ", name, "\n\n")

    res_contrast[[name]] <-
      fgczgseaora::runGSEA(
        data = filtered_dd,
        fpath = name,
        ID_col = ID_col,
        score_col = fc_col,
        organism = organism,
        target = target,
        nperm = nperm,
        outdir = result_dir,
        interestGeneType = idtype,
        contrast_name = filtered_dd$contrast[[1]]
      )
  }
  res[[target]] <- res_contrast
}

print(summary(warnings()))

saveRDS(res, file.path(result_dir, "GSEA_Results.Rda"))
copy_gsea_report(result_dir)
rmarkdown::render(file.path(result_dir, "GSEA_Results_Overview.Rmd"),
                  params=list(GSEA = res),
                  output_file = "index.html",
                  output_dir = result_dir)



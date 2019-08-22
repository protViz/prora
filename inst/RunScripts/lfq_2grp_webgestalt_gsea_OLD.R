options(warn = -1)


suppressMessages(library(WebGestaltR))

suppressMessages(library(tidyverse))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(sigora))
library(GO.db)
library(slam)
library(fgczgseaora)
library(readr)


# Files -------------------------------------------------------------------

args = commandArgs(trailingOnly=TRUE)
print(args)


if(length(args) < 2) {
  stop("At least two argument must be supplied: grp2file and output_directory", call.=FALSE)
}
#args <- c("data/MQ-report-mouseLFQ_8wPostSCI_vs_Nv.txt", "gsea_results")
grp2report <- args[1]
result_dir <- args[2]

organism <- "hsapiens"
if(length(args) == 3){
  organism <- args[3]
}
organisms <- listOrganism(hostName = "http://www.webgestalt.org/", cache = NULL)

if(! organism %in% organisms){
  stop("organism : ",organism, "is not in the list of available organisms", paste(organisms, collapse=" ") )
}


# number of permutations
nperm <- 10
if (length(args)==4) {
  nperm <- args[4]
}

ID_col <- "TopProteinName"
fc_col <- "log2FC"

if (length(args) == 6) {
  # default output file
  ID_col <- args[5]
  fc_col <- args[6]
}

#######################################

cat("grp2report", grp2report , "\n",
    "result_dir",result_dir , "\n",
    "organism", organism, "\n",
    "nperm",nperm , "\n",
    "ID_col" , ID_col,"\n",
    "fc_col" , fc_col , "\n")


target_GSEA <- c(
  "geneontology_Biological_Process",
  "geneontology_Cellular_Component",
  "geneontology_Molecular_Function"
)

#target_SIGORA <- target_SIGORA[1]
# Parameters --------------------------------------------------------------

#
fpath_se <- tools::file_path_sans_ext(basename(grp2report))
odir <- file.path( result_dir , make.names(fpath_se))

if(!dir.exists(result_dir)){
  dir.create(result_dir)
}

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



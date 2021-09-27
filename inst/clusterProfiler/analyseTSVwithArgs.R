# Author W Wolski
# this script analyses all projects in file output_MaxQuant_info_WSpecies.csv

library(tidyverse)
library(readr)

args = commandArgs(trailingOnly = TRUE)

parameters <- list()
if (length(args) == 4) {
  parameters$clustalg <- args[2]
  parameters$outfolder <- args[1]
  parameters$peptide <- as.logical(args[3])
  parameters$JK <- as.logical(args[4])
  print(parameters)
}

tmp <- readr::read_csv("output_MaxQuant_info_WSpecies.csv")
tmp$project_ID <-  gsub("/srv/www/htdocs/(p[0-9]{1,5})/bfabric/.*","\\1",tmp$path_to_zip)

human <- tmp %>% filter(species == "human" & n_inputs >= 4)

if (Sys.info()["sysname"] == "Windows") {
  human$windowpaths <- gsub("/srv/www/htdocs", "y:", human$path_to_zip)
  RscriptExe <- "c:/Program Files/R/R-4.1.1/bin/Rscript.exe"
  stopifnot(file.exists(RscriptExe))
  rscriptlocation <- "runscript.R"
}else{
  human$windowpaths <- human$path_to_zip
  RscriptExe <- "/usr/bin/Rscript"
  rscriptlocation <- "/scratch/CLUSTERPROFILER/clusterprofileron33/runscriptLocal.R"
}


for (i in 1:nrow(human)) {
  path <- human$windowpaths[i]
  workunitid <- human$workunitid[i]
  projectid <- human$project_ID[i]
  if (!file.exists(path)) {
    print(path)
  } else if (TRUE) {
    arguments <- c(rscriptlocation,
                   path,
                   "human",
                   parameters$outfolder,
                   parameters$clustalg,
                   workunitid,
                   projectid,
                   parameters$peptide,
                   parameters$JK)

      res <- system2(RscriptExe, args = arguments)
      if (res > 0) {
		print("ERROR START running script failed")
		cat(paste("ERROR : ", arguments ,"\n"))
		print("ERROR END")
	}
  }
}

#path <- "y:/p2621/bfabric/Proteomics/MaxQuant/2018/2018-10/2018-10-05/workunit_175797/721705.zip"
#path <- "y:/p2673/bfabric/Proteomics/MaxQuant/2021/2021-01/2021-01-19/workunit_256172/1792120.zip"
#path <-  "y:/p2954/bfabric/Proteomics/MaxQuant/2019/2019-01/2019-01-09/workunit_187912/1068470.zip"
#path <- "y:/p2865/bfabric/Proteomics/MaxQuant/2020/2020-01/2020-01-13/workunit_234979/1462850.zip"
if(FALSE){
path <- "/srv/www/htdocs/p2702/bfabric/Proteomics/MaxQuant/2019/2019-03/2019-03-14/workunit_194381/1163583.zip"
#path <- "y:/p2961/bfabric/Proteomics/MaxQuant/2019/2019-06/2019-06-23/workunit_200319/1292787.zip"
arguments <- c("runscriptLocal.R", path, "human","test" ,"DPA", "194381", "p2702", TRUE, FALSE)

commandArgs <- function(...){ return(arguments[-1]) }
commandArgs()
source("runscriptLocal.R")
}


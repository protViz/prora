library(tidyverse)
library(readr)
tmp <- readr::read_csv("output_MaxQuant_info_WSpecies.csv")
tmp$project_ID <-  gsub("/srv/www/htdocs/(p[0-9]{1,5})/bfabric/.*","\\1",tmp$path_to_zip)

human <- tmp %>% filter(species == "human" & n_inputs > 4)
summary(human$n_inputs)
hist(human$n_inputs)

human$windowpaths <- gsub("/srv/www/htdocs","y:",human$path_to_zip)

for (i in 1:nrow(human)) {
  path <- human$windowpaths[i]
  workunitid <- human$workunitid[i]
  projectid <- human$project_ID[i]
  if (!file.exists(path)) {
    print(path)
  } else if (TRUE) {
    arguments <- c("runscript.R", path, "human","blubA" ,"hclustdeepsplit", workunitid, projectid)
    res <- system2("c:/Program Files/R/R-4.0.3/bin/Rscript.exe", args = arguments)
    stopifnot(res == 0)
    if (TRUE) {
      commandArgs <- function(...){return(arguments[-1])}
      commandArgs()
      source("runscript.R")
    }
  }
}

path <- "y:/p2621/bfabric/Proteomics/MaxQuant/2018/2018-10/2018-10-05/workunit_175797/721705.zip"
path <- "y:/p2673/bfabric/Proteomics/MaxQuant/2021/2021-01/2021-01-19/workunit_256172/1792120.zip"
path <-  "y:/p2954/bfabric/Proteomics/MaxQuant/2019/2019-01/2019-01-09/workunit_187912/1068470.zip"
path <- "y:/p2865/bfabric/Proteomics/MaxQuant/2020/2020-01/2020-01-13/workunit_234979/1462850.zip"
arguments <- c("runscript.R", path, "human","blubB" ,"hclustdeepsplit", "p3220", "23333")

commandArgs <- function(...){return(arguments[-1])}
commandArgs()
source("runscript.R")



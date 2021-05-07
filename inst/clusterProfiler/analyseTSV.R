library(tidyverse)
library(readr)
tmp <- readr::read_csv("output_MaxQuant_info_WSpecies.csv")
human <- tmp %>% filter(species == "human" & n_inputs > 4)

summary(human$n_inputs)
hist(human$n_inputs)
windowpaths <- gsub("/srv/www/htdocs","y:",human$path_to_zip)

for (path in windowpaths[1:length(windowpaths)]) {
  if (!file.exists(path)) {
    print(path)
  } else if (TRUE) {
    arguments <- c("runscript.R", path, "human","blubC" ,"hclustdeepsplit")
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
arguments <- c("runscript.R", path, "human","blubB" ,"hclustdeepsplit")

commandArgs <- function(...){return(arguments[-1])}
commandArgs()
source("runscript.R")



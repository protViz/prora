library(tidyverse)
library(readr)
tmp <- readr::read_csv("output_MaxQuant_info_WSpecies.csv")
names(tmp)
human <- tmp %>% filter(species=="human" & n_inputs > 4)

summary(human$n_inputs)
hist(human$n_inputs)
windowpaths <- gsub("/srv/www/htdocs","y:",human$path_to_zip)

for (path in windowpaths[4:length(windowpaths)]) {
  if (!file.exists(path)) {
    print(path)
  }else if(TRUE){
    arguments <- c("runscript.R", path, "human","blubA" ,"hclust")
    res <- system2("c:/Program Files/R/R-4.0.3/bin/Rscript.exe", args = arguments)
    stopifnot(res == 0)
    if (FALSE) {
      commandArgs <- function(...){return(arguments[-1])}
      commandArgs()
      source("runscript.R")
    }
  }
}



library(tidyverse)
library(readr)

path <- "20200511_hclust_V1"
allFiles <- dir(path)
summaryFiles <- grep("Summary_", allFiles, value = TRUE)


res <- vector(mode = "list", length = length(summaryFiles))

for (i in 1:length(xx)) {
  print(summaryFiles[[i]])

  res[[i]] <- read_tsv(file = file.path(path,summaryFiles[i]))
}

summary <- bind_rows(res)
head(summary)
table(summary$id.mapping.service)

pairs(summary[,7:(ncol(summary)-1)])
s3500E <- summary %>% filter(nr.ENTREZIDS > 3500  & nr.ENTREZIDS < 4000)
pairs(s3500E[,7:(ncol(summary)-1)])
s3500E <- summary %>% filter(nr.ENTREZIDS > 3500  & nr.ENTREZIDS < 4000)

View(s3500E)

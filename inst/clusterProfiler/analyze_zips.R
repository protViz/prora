library(readr)
library(tidyverse)
tmp <- readr::read_csv("output_MaxQuant_info.csv")

sort(table(tmp$fastaFiles))


tmp2 <- tmp %>% mutate(species = case_when(grepl("human",fastaFiles) ~ "human",
                                   grepl("9606", fastaFiles) ~ "human",
                                   grepl("10090", fastaFiles) ~ "mouse",
                                   grepl("4932", fastaFiles) ~ "yeast",
                                   grepl("9913", fastaFiles) ~ "bowine",
                                   grepl("10116", fastaFiles) ~ "rat",
                                   grepl("9615", fastaFiles) ~ "dog",
                                   TRUE ~ fastaFiles))
sort(table(tmp2$species))

write_csv(tmp2, file="output_MaxQuant_info_WSpecies.csv")

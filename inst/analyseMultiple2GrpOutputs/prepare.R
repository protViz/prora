rm(list = ls())
library(tidyverse)
library(prora)
library(msigdbr)

path <- "data"

GSIDX  <- 1
outdir <- "newdelivery"



files <- dir(path, pattern = "*.txt")
res <- vector(mode = "list", length = length(files))
for (i in 1:length(files)) {
  file <- files[[i]]
  data <- read_tsv(file.path(path, file))
  data$contrast <- file
  data <- data %>% select(contrast ,TopProteinName, log2FC, pseudo.log2FC, t)
  res[[i]] <- data
}
res <- bind_rows(res)
res <- res %>% mutate(contrast = gsub(".txt","", contrast))
res <- prora::get_UniprotID_from_fasta_header(res, idcolumn = "TopProteinName")
res <- prora::map_ids_uniprot(res)

write_tsv(res, file = file.path(path,"allContrasts.txt"))


res <- na.omit(res)
# summarize statistics for multiple mappings
res <- res %>% group_by(contrast, P_ENTREZGENEID) %>% summarize(log2FC = mean(log2FC), pseudo.log2FC = mean(pseudo.log2FC), t = mean(t))
ranklist <- fgsea_rank_contrasts(res, ids = "P_ENTREZGENEID", score = "pseudo.log2FC"  )



# Prepare gene sets

msigdbr::msigdbr_species()
species <- "Homo sapiens"
C5 <- msigdbr_collections() %>% filter(gs_cat == "C5")
C5 <- msigdbr_collections() %>% filter(gs_cat == "H")

getMsigdbGenesets <- function(msigCollection, species){
  genesetsC5 <- vector(mode = "list", length = nrow(msigCollection) )
  for (i in 1:nrow(msigCollection)) {
    genesetsC5[[i]] <- msigdbr(species = species,
                               category = msigCollection$gs_cat[i],
                               subcategory = msigCollection$gs_subcat[i])
  }
  names(genesetsC5) <- make.names(msigCollection$gs_subcat)
  fgseaGSlist <- lapply(genesetsC5 , fgsea_msigdb)
  return(fgseaGSlist)
}
fgseaGSlist <- getMsigdbGenesets(C5, species)



fgseaRes <- prora::run_fgsea_for_allContrasts(ranklist , fgseaGSlist[[GSIDX]])



all <- bind_rows(fgseaRes)
all <- all %>%
  dplyr::relocate(nMoreExtreme, pval,  ES, leadingEdge, .after = size) %>%
  dplyr::relocate(comparison ,.before = pathway)

fgseaResultsForAllContrasts <- paste0("Geneset_",names(fgseaGSlist)[GSIDX],"_allContrasts.xlsx")
all <- tibble::as_tibble(all)
all$leadingEdge = sapply(all$leadingEdge, function(x){paste(x, collapse = ";")})
writexl::write_xlsx(all, path = file.path(outdir, fgseaResultsForAllContrasts))

markdown <- paste0("Geneset_",names(fgseaGSlist)[GSIDX],"_allContrasts.html")
rmarkdown::render("multigroupGSEA.Rmd", params = list(all = all, threshold = 0.2), output_format = bookdown::html_document2())
file.copy("multigroupGSEA.html" , file.path(outdir, markdown), overwrite = TRUE )

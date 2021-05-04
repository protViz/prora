rm(list = ls())
library(tidyverse)
library(prora)
library(msigdbr)

path <- "data2grp/"

threshold = 0.2

GSIDX  <- 1
outdir <- "newdeliveryOther"
dir.create(outdir)


files <- dir(path, pattern = "*.txt")
res <- vector(mode = "list", length = length(files))

for (i in 1:length(files)) {
  file <- files[[i]]
  data <- read_tsv(file.path(path, file))
  data$contrast <- file
  data <- data %>% select(contrast , TopProteinName, log2FC, pseudo.log2FC, t)
  res[[i]] <- data
}

res <- bind_rows(res)
res <- res %>% mutate(contrast = gsub(".txt","", contrast))
res <- prora::get_UniprotID_from_fasta_header(res, idcolumn = "TopProteinName")
res <- prora::map_ids_uniprot(res)
res <- res %>% filter(grepl("^MQ_*", contrast))
write_tsv(res, file = file.path(outdir , "allContrasts.txt"))


res <- na.omit(res)
# summarize statistics for multiple mappings
res <- res %>% group_by(.data$contrast, .data$P_ENTREZGENEID) %>%
  summarize(log2FC = mean(.data$log2FC), pseudo.log2FC = mean(.data$pseudo.log2FC), t = mean(.data$t))

ranklist <- prora::fgsea_rank_contrasts(res, ids = "P_ENTREZGENEID", score = "t"  )


# Prepare gene sets
msigdbr::msigdbr_species()
species <- "Homo sapiens"
species <- "Mus musculus"

hallmark <- {msigdbr_collections() %>% filter(.data$gs_cat == "H")}

#hallmark$gs_subcat <- "HALLMARK"
C5 <- bind_rows( {msigdbr_collections() %>% filter(.data$gs_cat == "C5") %>% filter(grepl("^GO:", .data$gs_subcat))},
                 hallmark,
                 {msigdbr_collections() %>% filter(.data$gs_subcat == "CP:KEGG")} )


fgseaGSlist <- prora::getMsigdbGenesets(C5, species)
fgseaRes <- prora::run_fgsea_for_allContrasts(ranklist , fgseaGSlist[[GSIDX]],minSize = 15)

all <- bind_rows(fgseaRes)
all <- all %>%
  dplyr::relocate(nMoreExtreme, pval,  ES, leadingEdge, .after = size) %>%
  dplyr::relocate(comparison ,.before = pathway)

fgseaResultsForAllContrasts <- paste0("Geneset_",names(fgseaGSlist)[GSIDX],"_allContrasts.xlsx")
all <- tibble::as_tibble(all)
all$leadingEdge = sapply(all$leadingEdge, function(x){paste(x, collapse = ";")})
writexl::write_xlsx(all, path = file.path(outdir, fgseaResultsForAllContrasts))


## Move logic from markdown

xxNES <- tidyr::pivot_wider(all, id_cols = "pathway", names_from = "comparison" , values_from  = "NES")
minpadj <- all %>% group_by(pathway) %>% summarise(minpadj = min(padj, na.rm = TRUE))
xxNES <- inner_join(xxNES, minpadj)
dim(xxNES)

xxNES <- xxNES %>%  filter(minpadj < threshold)  %>% arrange(minpadj)

xmNES <- xxNES %>% dplyr::select( -all_of(c("pathway", "minpadj") )) %>% as.matrix
rownames(xmNES) <- xxNES$pathway


### End of


markdown <- paste0("Geneset_",names(fgseaGSlist)[GSIDX],"_allContrasts.html")
rmarkdown::render("multigroupGSEA.Rmd",
                  params = list(allpvalues = all,
                                minpvalue = minpadj,
                                xmNES = xmNES ,
                                xxNES = xxNES, gsName = names(fgseaGSlist)[GSIDX]),
                  output_format = bookdown::html_document2())

file.copy("multigroupGSEA.html" , file.path(outdir, markdown), overwrite = TRUE )


